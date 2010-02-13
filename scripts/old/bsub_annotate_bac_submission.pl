#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Temp qw/tempdir/;
use File::Basename;

use Data::Dumper;
use CXGN::Publish;
use CXGN::TomatoGenome::BACSubmission;
use CXGN::TomatoGenome::BACPublish qw/publishing_locations bac_publish resource_file cached_validation_text/;
use CXGN::Tools::Wget qw/wget_filter/;
use CXGN::IndexedLog;

########### DEFAULTS ##########

my $conf = CXGN::TomatoGenome::Config->load_locked;
my $default_publish_dir = File::Spec->catdir($conf->{'ftpsite_root'},$conf->{'bac_publish_subdir'});
my $default_temp_dir    = File::Spec->tmpdir;
my $default_log = $conf->{'bac_processing_log'};

###############################

#some debugging routines
use constant DEBUG => ($ENV{BSUBANNOTDEBUG} ? 1 : 0);
sub dprint(@) { if(DEBUG) { local $|=1; print STDERR @_; } }
sub vprint(@) { print STDERR @_ if DEBUG || our $verbose; }

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] <bac submission tarball> ...

  Run a suite of automatic annotation tools on the given BAC tarball(s),
  publishing the annotation results in either the same directory as
  the tarball (the default)

  Options:

  -d <dir>
       The filesystem directory to publish the result files to.
       Default: $default_publish_dir

  -t <dir>
       Directory for storing temp files.  Usually needs 500-700MB.
       Default: $default_temp_dir

  -f ignore validation status of the given submission.  Try to annotate
     regardless.

  -l <table>
     Log processing messages to given log table.
     Default: $default_log

  -x Dry run.  Don't actually copy or publish anything.

  -v Be verbose (to STDERR)

EOU
}

###get and validate our command line options
my %opt;
getopts('d:t:l:xvf',\%opt) or usage();
@ARGV or usage('no submission tarballs given');
#d
$opt{d} ||= $default_publish_dir;
-d $opt{d} or die("publish dir '$opt{d}' does not exist");
unless(-w $opt{d} || $opt{x}) {
  system "ls $opt{d}";
  -w $opt{d} or die("publish dir '$opt{d}' is not writable");
}
#t
if( $opt{t} ) {
    $opt{t} = File::Spec->rel2abs($opt{t});
    CXGN::Tools::Wget->temp_root_dir($opt{t});
}
$opt{t} ||= $default_temp_dir;
-d $opt{t} or die("temp dir '$opt{t}' does not exist");
-w $opt{t} or die("temp dir '$opt{t}' is not writable");
CXGN::TomatoGenome::BACSubmission->tempdir($opt{t});
#l
$opt{l} ||= $default_log;
my $dbh = CXGN::DB::Connection->new;
my $log = CXGN::IndexedLog->open(DB => $dbh, $opt{l} );
#x
$CXGN::Publish::dry_run = 1 if $opt{x};
#v
if( $opt{v} ) {
  our $verbose = 1;
  $CXGN::Publish::print_ops = 1;
}

vprint "fetching resource files...\n";
### set up some parameters for the analyses we'll be running
our $anal_params = {
		    repeatmasker_lib_file  => resource_file('repeats_master'),
		    geneseqer_est_seq_file => resource_file('sgn_ests_tomato'),
		    geneseqer_ug_seq_file  => resource_file('lycopersicum_combined_unigene_seqs'),
		    gth_sgne_seq_file      => resource_file('sgn_ests_tomato'),
		    gth_sgnu_seq_file      => resource_file('lycopersicum_combined_unigene_seqs'),
		   };
vprint "done fetching resource files.\n";


### loop over each submission and annotate it
my @failed_files;
my @successful_files;
foreach my $file (@ARGV) {
  my $file_done = 0;
  eval {
    vprint "beginning annotation for $file...\n";
    $log->append("starting annotation for ".basename($file)) if $log;
    annotate_submission($file);
    vprint "finished annotating $file.\n";
    $file_done = 1;
  }; if(!$EVAL_ERROR && $file_done) {
    push @successful_files,$file;
  } else {
    my ($basename) = fileparse($file);
    my $errstr = "$basename processing failed:\n$EVAL_ERROR\n";
    warn $errstr;
    $errstr =~ s/\n/\\n/g; #replace any newlines
    $log->append($errstr) if $log;
    push @failed_files,$file;
  }
  print "Finished $file.\n";
}
print "Succeeded: ",join(', ',@successful_files),"\n";
print "Failed: ",join(', ',@failed_files),"\n";


#here's where the real work gets done.  open the submission,
#run the annotation pipeline on it, then publish the results
sub annotate_submission {
  my $filename = shift;

  #before opening it, check whether it's already known to be invalid
  vprint "checking for cached validation errors...\n";
  my $cached_validation_text = cached_validation_text($filename);
  if(!$opt{f} and $cached_validation_text) {
    die "from validation cache: $cached_validation_text";
  }
  vprint "done checking for cached validation errors.\n";

  vprint "Opening $filename...\n";

  my $sub = CXGN::TomatoGenome::BACSubmission->open_stripped($filename);
  vprint "done.\n";

  #only validate the submission if it's not in the cache
  #the validation process will update the cache
  if(!$opt{f} and !defined($cached_validation_text) and my $errortext = $sub->validation_text ) {
    die $errortext;
  }

  #now run the annotation pipeline on it
  vprint "Submission validated, analyzing submission...\n";
  my %analysis_result_files = $sub->analyze_new_submission($anal_params);
  vprint "done.  Got result files:\n",Dumper(\%analysis_result_files);

  #figure out the publishing operations needed to publish the result files
  my @publish_operations;
  #$destinations is a hashref of full paths for publishing the result files
  my $destinations = publishing_locations($opt{d},$sub->sequence_identifier,$sub->is_finished);
  #for each analysis,
  while ( my( $analysis_name, $analysis_results ) = each  %analysis_result_files ) {
    #if it has results
    if ( $analysis_results && @$analysis_results ) {
      #look up the destinations for this analysis
      my $dest = $destinations->{"annot_$analysis_name"}
	or die "No destination set for result files from analysis '$analysis_name'";
      #and anywhere old versions of these might be published
      my $old_dest = $destinations->{obsolete}->{"annot_$analysis_name"};

      #for each result file, look up its destination, the publishing
      #destination for any old version of it, and make two publish operations,
      #one to publish this result file, and another to unpublish its
      #old version (in the unfinished directory, if this BAC is finished,
      #             or in the finished directory, if this BAC is unfinished)
      #if the unfinished version isn't there, the publish will just skip removing it.
      #that's what the 'rm -f' publishing opcode means

      #check the lengths of the arrays and snark if they're not right
      @$analysis_results == @$dest
	or die "Must have as many publishing destinations as we have result files.\n",
	  "I have results:\n",Dumper($analysis_results),"and I have copy destinations:\n",
	    Dumper($dest),"and they don't match up.\n";
      @$analysis_results <= @$old_dest
	or die "Must have at least as many alternate publishing destinations as we have result files.\n",
	  "I have results:\n",Dumper($analysis_results),"and I have old version publish locations:\n",
	    Dumper($old_dest),"and they don't match up.\n";

      #now make the proper publish operations
      for(my $i=0; $i<@$analysis_results; $i++) {
	push @publish_operations,[ 'cp',    $analysis_results->[$i],$dest->[$i] ];
      }
      push @publish_operations, map {['rm -f', $_]} @$old_dest;
    }
  }

  vprint "Assembled publish operations:\n",Dumper(\@publish_operations);

  #now do the actual publish, which will die if it fails
  local $CXGN::Publish::make_dirs = 1; #create any destination dirs that don't exist
  bac_publish(@publish_operations);

  #write our success to the log file if one was specified
  $log->append("ANNOTATION_COMPLETED_FOR ".basename($filename)) if $log;
}


