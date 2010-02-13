#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Basename;

use CXGN::DB::Connection;

use CXGN::Tools::Run;

use CXGN::TomatoGenome::BACPublish qw/publishing_locations parse_filename/;
use CXGN::Publish qw/mkdir_or_print/;
use CXGN::IndexedLog;

use CXGN::TomatoGenome::Config;
use Hash::Util qw/lock_hash/;

use Data::Dumper;


######## DEFAULTS ###########

my $cfg = CXGN::TomatoGenome::Config->load_locked;
lock_hash(%$cfg);
my $default_log = $cfg->{'bac_processing_log'};
my $default_job_logs_dir = File::Spec->catdir($cfg->{'bac_pipeline_dir'}, $cfg->{'bac_job_logs_dir'});

#############################

#note that these debugging routines return their arguments
use constant DEBUG => $ENV{BSUBANNOTATEALLDEBUG} ? 1 : 0;
$ENV{CXGN_DEBUG} = 1 if DEBUG;
sub vprint(@) {
  print @_ if our $verbose or our $dry_run;
  return @_;
}
sub dprint(@) {
  print @_ if DEBUG;
  return @_;
}
sub dprinta(@) {
  dprint join(' ',@_),"\n";
  return @_;
}


sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] <publishing directory>
  $FindBin::Script [options] -r

  Script to look in the given publishing directory, find BAC submissions
  that need to be annotated or re-annotated, and submit jobs to the CXGN
  cluster for annotating them.  Note that cluster nodes must be able to
  write directly to this directory.

  Arguments:

  <publishing directory>
      Publishing directory to search for annotatable BACs.  Worker nodes
      must be able to write directly to this directory.

  -r  Reset the annotation pipeline and exit.  The next run of the pipeline
      will result in annotations being rerun on ALL BACs in the publishing
      directory.

  -l <table>
      Database table to use for writing processing logs, which are also used
      as the persistent state of the annotation pipeline.  Table must be
      writable.
      Default: $default_log

  -L  <dir>
      Directory for per-job logs, etc.
      Default: $default_job_logs_dir

  -f  Force attempted annotation of invalid (by current standards)
      submissions

  -x  Dry run.  Don't submit any jobs or write any persistent state, just
      print what operations would be performed.  Implies -v.

  -q  <queue>
      submit the jobs to the specified torque queue.
      defaults to 'batch\@solanine.sgn.cornell.edu'

  -v  Be verbose;

EOU
}

our %opt;
getopts('q:l:L:rxvf',\%opt) or usage();
$opt{L} ||= $default_job_logs_dir;
$opt{q} ||= 'batch@solanine.sgn.cornell.edu';
if($opt{x}) { #set dry run flags
  our $dry_run = 1;
  $CXGN::Publish::dry_run = 1;
}
if($opt{v}) { #set verbose flags
  our $verbose = 1;
  $CXGN::Publish::print_ops = 1;
}

#check that we have the necessary binaries in our path
foreach my $bin (qw/find qsub qstat/) {
  `which $bin` or die "Could not find '$bin' in your path.  Is it installed?\n";
}

#$pubdir is the path to the publication directory (e.g. ftp site)
my $pubdir = shift @ARGV;
@ARGV and usage();
unless($opt{r}) {
  -d $pubdir or $opt{r} or  usage("no directory '$pubdir' found");
  -w $pubdir or $opt{r} or  usage("'$pubdir' not writable");
}

#$log_dir is the path to the directory where we keep all our
#annotation logs and pipeline state
my $annotation_job_logs_dir = $opt{L};
-d $annotation_job_logs_dir
  or mkdir_or_print($annotation_job_logs_dir)
  or die "Could not create job logs directory '$annotation_job_logs_dir': $!";

#$logtable_name is the name of the annotation log table,
#which keeps state for the annotation pipeline
my $dbh = CXGN::DB::Connection->new;
my $logtable_name = 'cxgn_bac_pipeline_processing_log';

#append a truncation mark to the log file if the -r option was passed
#don't ever _really_ truncate it, since we might want to be able to easily
#undo a pipeline reset
if($opt{r}) {
  print "Resetting annotation pipeline...";
  unless(our $dry_run) {
    my $log = CXGN::IndexedLog->open(DB => $dbh, $logtable_name);
    $log->reset;
    print "done.\n";
  } else {
    print "dry run, not actually done.\n";
  }
  exit;
}

### look in the publishing directory, find all tar files in the chr## dirs
### with well-formed names.
my @chrdirs = (map { "$pubdir/".sprintf("chr%02d",$_) } (0..12));
vprint "Searching publication dir for BAC submissions...\n";
my @submissions = map {
  chomp;
  if(my $p = parse_filename($_)) {
    ($p)
  } else {
    warn "Tarfile $_ does not appear to be a BAC submission.  Skipping.\n";
    ()
  }
} map glob("$_/*/*.tar.gz"), @chrdirs;
vprint "done.\n";

dprint "got submissions:\n",Dumper(\@submissions);

@submissions or die "No BAC submission files found in '$pubdir'";

# now read in the annotation file, storing in %annot_status hash,
# but clear the annot_status every time we hit a RESET,
# so that at the end the status only has whatever was
# after the last reset.
vprint "Reading log in '$logtable_name'...";
my $log = CXGN::IndexedLog->open(DB => $dbh, $logtable_name);
$log->lookup(content => 'foobar baz'); #< do a lookup to force it to parse the log
vprint "done.\n";

### now submit annotation jobs for all the ones that need it

foreach my $submission (@submissions) {
  if(my %submissionrecord  = $log->lookup( content => "ANNOTATION_COMPLETED_FOR $submission->{basename}")) {
    vprint "$submission->{basename} - up to date\n";
    next;
  } elsif( %submissionrecord = $log->lookup( content => "SUBMITTED_ANNOTATION_JOB_FOR $submission->{basename}") ) {
    #check if the job is still in the queue.  if so, skip.  if not, it must have died or something,
    #so we'll continue on and submit another one
    my (undef,undef,$jobid) = split /\s+/,$submissionrecord{content};
    $jobid =~ s/(\\n)+$//;
    $jobid =~ /^\d+(\.\w+)*$/ or die "Invalid jobid '$jobid' found in log file";
    if(`qstat $jobid 2>/dev/null`) {
      vprint "$submission->{basename} - already has an annotation job in the queue ($jobid), not queueing another.\n";
      next;
    } else {
      vprint "$submission->{basename} - $jobid seems to have died.  submitting another.\n";
    }
  }

  vprint "$submission->{basename} - submitting annotation job\n";
#   my $pubfiles = publishing_locations($pubdir,$submission->{seq_name},$submission->{finished});
#   #take the annotation files in the pubfiles hash, flatten them into one long array
#   my @annotfiles =
#     map  { s|^$pubdir/*||; $_} #remove the publication directory from them
#     map  { @$_ } #flatten them into one array, since they are each arrayrefs
#     @{$pubfiles}{ grep {/^annot_/} keys %$pubfiles }; #get all the pubfile entries that have keys beginning in 'annot_'

  #submit a job to the cluster to annotate this BAC
  my $f = $opt{f} ? '-f' : ''; #< try to annotate invalid submissions if this is given
  my $annot_cmd = "bsub_annotate_bac_submission.pl -d $pubdir $f -l $logtable_name $submission->{filename}";
  dprint "annotation command: '$annot_cmd'\n";

  my $jobid = '0.fakejob.id';
  #this is where the individual annotation jobs will put their output and such

  unless( our $dry_run ) {
    CXGN::Tools::Run->run(
			  dprinta
			  qsub => ( -N => "an$submission->{seq_name}",
				    -j => 'oe', #join the stdout and stderr of the job
				    -o => $annotation_job_logs_dir, #keep the job's result logs
				    -d => '/', #set its working dir to /, to avoid any confusion about paths not being symlinked the same on different machines
				    -q => $opt{q},
                                    -l => 'nodes=1,vmem=2000M',
				  ),
			  { in_file  => \$annot_cmd,
			    out_file => \$jobid,
			  }
			 )
	;
    $log->append(dprinta("SUBMITTED_ANNOTATION_JOB_FOR $submission->{basename} $jobid"));
    vprint "Submitted annotation job $jobid\n";
    sleep 1; #don't submit jobs too fast, it confuses the scheduler
  } else {
    vprint "Skipping submission of annotation job (dry run)\n";
  }
}


