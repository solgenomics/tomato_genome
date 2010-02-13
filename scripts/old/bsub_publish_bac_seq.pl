#!/usr/bin/perl
use strict;
use warnings;
use English;
use Carp;
use LWP::Simple; #for downloading resource files

use FindBin;

use File::Basename;
use File::Copy;
use File::Spec;
use File::Find;
use File::Path;
use File::Temp qw/tempfile tempdir/;

use Data::Dumper;

use Memoize;

use Bio::Index::Fasta;

use CXGN::Tools::File qw/ file_contents size_changing /;
use CXGN::Publish qw/publish published_as move_or_print orig_basename/;

use Getopt::Std;
use CXGN::TomatoGenome::BACSubmission qw/:errors/;
use CXGN::TomatoGenome::BACPublish qw/ parse_filename publishing_locations bac_publish aggregate_filename cached_validation_text cached_validation_errors /;
use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident assemble_clone_ident/;
use CXGN::IndexedLog;
use CXGN::Tools::Script qw/lock_script unlock_script/;

############### CONFIGURATION VARIABLES ##################

use CXGN::TomatoGenome::Config;
my $conf = CXGN::TomatoGenome::Config->load_locked;
my $country_uploads_path = $conf->{'country_uploads_path'};
my $publish_path = File::Spec->catdir($conf->{'ftpsite_root'},$conf->{'bac_publish_subdir'});
my $default_log = $conf->{'bac_processing_log'};

##########################################################

#some debugging routines
use constant DEBUG => ($ENV{BSUBPUBLISHDEBUG} ? 1 : 0);
sub dprint(@) { if(DEBUG) { local $|=1; print STDERR @_; } }

#make a lock file to keep multiple instances of this program from running
lock_script() or die "Don't run more than one instance of $FindBin::Script";

sub usage {
  unlock_script();
  die <<EOU;
  Usage:
    $FindBin::Script [-d <dirname>] [-rx] <submission tarballs>

  Script to take a BAC submission tarball, vector screen its sequence file,
  and copy it to the publishing location (e.g. the FTP site directories).
  If the submission tarball(s) are in a country upload directory
  ($country_uploads_path/<country>/upload), the script will attempt,
  after processing, to move the tarball to upload-processed/ in that country's
  home dir.


  Options:
   -d <dirname>
        publishing base directory
        (default $publish_path)
   -s <sp_organization_id>
        force script to set sequencer info using the given sp_organization_id,
        from the sgn_people.sp_organization.sp_organization_id column in the DB
   -r   attempt to repair improperly formatted bac submission tarballs
   -c   enable caching of validation results
   -x   dry run, do not modify or move anything, just print what you would do.
   -l <table name>
        log operations to the given table in the cxgn database.
        Defaults to 'bac_processing_log' conf variable, which is currently set
        to '$default_log'
   -f   fast, don't wait 70 seconds to see if file sizes are changing.
   -S   after publishing everything given on the command line, go through the
        BAC repository and update file names and contents to make sure all the
        CXX* file names agree with the chromosome assignments in the BAC
        registry

   -G   accept submissions with missing or malformed GenBank accessions or entries
   -M   accept submissions with multiple sequences (disallowed by current guidelines)
   -V   accept submissions with vector sequence content

   NOTE: to change the database host, user, and password, set environment
         variables DBHOST, DBUSER, and DBPASS, or change your configuration in
         /etc/cxgn/Default.HostConf.
EOU
}

our %opt; #GLOBAL options
getopts('d:s:crfxl:GMSV',\%opt) or usage;
$CXGN::Publish::dry_run = 1 if $opt{x};
$CXGN::Publish::print_ops = 1 if DEBUG;

my @submission_files = @ARGV;
@submission_files || $opt{S} or usage;

$opt{d} ||= $publish_path;
($opt{d} && -d $opt{d}) || $opt{x} or die "Destination dir '$opt{d}' not found\n";
($opt{d} && -w $opt{d}) || $opt{x} or die "Destination dir '$opt{d}' not writable\n";
$opt{l} ||= $default_log;

my $dbh = CXGN::DB::Connection->new;
my $log = CXGN::IndexedLog->open(DB => $dbh,  $opt{l});
$log->is_writable or die "log '$opt{l}' not writable!";

my @failed_files;
my @successful_files;

dprint "idents index is ",Dumper _identifiers_index() if DEBUG;

#for each submission file we've been given
foreach my $file (@submission_files) {

  #NOTE: process_file raises an error if it fails
  print "Processing $file...\n";
  my $file_done = 0;
  eval {
    $log->append("starting publish for ".basename($file)) if $log;
      #check the file to make sure its size is not changing right now
      #which would mean it is still uploading
      size_changing($file,$opt{f} ? 70 : 1)
	  and die "Skipping '$file', its size is still changing.  Still being uploaded?\n";
      process_file($file);
      $file_done = 1;
  }; if(!$EVAL_ERROR && $file_done) {
    push @successful_files,$file;
    $log->append("successfully published $file") if $log;
  } else {
    my ($basename) = fileparse($file);
    my $errstr = "$basename processing failed:\n$EVAL_ERROR\n";
    print $errstr;
    $errstr =~ s/\n/\\n/g;
    $log->append($errstr) if $log;
    push @failed_files,$file;
  }
  print "Finished $file.\n";
}

print "Submissions that passed:\n", map {"$_\n"} @successful_files;
print "Submissions that failed:\n", map {"$_\n"} @failed_files;

#if the -S option was passed, run the procedure to sync BAC names to chromosome nums
if( $opt{S} ) {
  sync_bac_names_in_repository(\@successful_files);
}

unlock_script(); #remove program lock file

#exit with error if any files failed
exit @failed_files ? 1 : 0;

##################### subroutines ##########################

#process a single submission file, which involves untarring it,
#looking at it, maybe repairing it, analyzing it in various ways,
#and finally copying the results of all these things to the
#publishing directory.
#if anything goes wrong, throws an error with die()

sub _unallowable_errors; #< predeclare _unallowable_errors sub

sub process_file {
  my $submission_file = shift;
  my ($submission_file_basename,$submission_file_path) =
    fileparse($submission_file);
  $submission_file_path = File::Spec->rel2abs($submission_file_path);

  #before opening, check if we already know it won't validate
  if( ! $opt{r} && $opt{c} && _unallowable_errors cached_validation_errors($submission_file) ) {
    die "from cache: ".cached_validation_text($submission_file);
  }

  #open the tar file
  -r $submission_file || die "Could not open submission file '$submission_file'\n";
  my $submission = CXGN::TomatoGenome::BACSubmission->open($submission_file);

  $submission->repair if $opt{r};

  ### die informatively if the submission has non-allowable errors
  my @errs = _unallowable_errors $submission->validation_errors;
  if( @errs ) {
    my $repairstr = $opt{r} ? "Repair failed. " : '';
    die $repairstr, $submission->validation_text;
  }

  #make sure we haven't already processed this bac once this session.
  #we can't do a given bac more than once per run, because then the
  #indexes of the sequence versions wouldn't work, and that would not be cool.
  our %bacs_already_processed;
  if($bacs_already_processed{$submission->bac_name}) {
    die "Already processed a BAC named ".$submission->bac_name." during this script run.  Skipping.";
  } else {
    $bacs_already_processed{$submission->bac_name} = 1;
  }

  #since this file is valid if we have gotten here,
  #figure out what upload account it comes from and set it in the tarball
  if( !$opt{s} &&  $submission_file_path =~ m|$country_uploads_path/([^/]+)/upload(-processed)?/*$| ) {
    -f $submission_file or die 'Sanity check failed';

    # check whether the submission already has sequencer information
    # in sequencer_info.txt and only rewrite it if not
    unless( $submission->sequencing_info ) {
      $submission->sequencing_info( org_upload_account_name => $1 );
    }
  }
  elsif( $opt{s} && !$submission->sequencing_info ) {
      $submission->sequencing_info( org_sp_organization_id => $opt{s} + 0);
  }

  #find and set the sequence version of this submission
  set_submission_sequence_version($submission);

  #generate file of renamed, vector-screened sequences
  my $seqs_file = $submission->vector_screened_sequences_file;

  #overwrite its old sequences file with the renamed sequences in temp storage
  copy($submission->renamed_sequences_file,$submission->sequences_file)
    or die "Could not copy within temp storage: cp ".$submission->renamed_sequences_file.' '.$submission->sequences_file.": $!\n";

  #tar up the newly relabeled and version-numbered and repaired stuff
  my $new_tarfile = $opt{x} ? 'new_tarball.tar.gz' : $submission->new_tarfile;

  {## last, if everything worked okay, move the information into the
   ## publishing directory (ftp site)

    my @publish_operations;
    #list of files where this submission should be published
    dprint "publishing submission with sequence identifier: ".$submission->sequence_identifier."\n";
    my $destinations = publishing_locations($opt{d},$submission->sequence_identifier,$submission->is_finished);

    #publish sequence file, and remove any obsoleted versions
    push @publish_operations,[ 'cp', $submission->vector_screened_sequences_file, $destinations->{seq} ];
    push @publish_operations, map {['rm -f',$_]} @{$destinations->{obsolete}{seq}};

    #publish (possibly repaired) tar file, and remove any otherly-finished version
    push @publish_operations,[ 'cp',    $new_tarfile, $destinations->{tar} ];
    push @publish_operations, map {['rm -f',$_]} @{$destinations->{obsolete}{tar}};

    dprint "doing publish operations:\n",Dumper(\@publish_operations);

    #now do the publish, which will die if it fails
    local $CXGN::Publish::make_dirs = 1; #make the directories needed during the publish
    bac_publish( @publish_operations );
  }

  #delete anything in this submission from temporary storage
  $submission->close;

  #if we've gotten here, everything went OK.  Move this file to its finished directory.
  if( $submission_file_path =~ m|$country_uploads_path/[^/]+/upload/*$| ) {
    $submission_file_path =~ s/\/upload(?=\/*$)/\/upload-processed/;
    print "Moving $submission_file_basename to $submission_file_path.\n";
    move_or_print( $submission_file, $submission_file_path)
      or die "Could not move $submission_file_basename from upload to upload-processed: $!\n";
  } else {
    print "Submission file $submission_file_basename is not in an upload directory, skipping move to upload-processed/.\n";
  }
}


sub sync_bac_names_in_repository {
  my ( $successful_submission_files ) = @_;

  #make a hash of clone names from the submission files
  my %passed_submissions = map {
    my $p = parse_filename($_);
    $p->{clone_name} => 1
  } @$successful_submission_files;

  my @tarfiles = glob File::Spec->catfile($opt{d},'chr*/*finished/*.tar.gz');

  # for each tarball in the repository
  foreach my $tarfile (@tarfiles) {

    # parse its filename
    my $p = parse_filename($tarfile);
    unless($p) {
      warn "WARNING: could not parse filename $tarfile";
      next;
    }

    # skip if there was a submission processed in this round for this
    # clone.  we don't want to overwrite whatever that had
    next if $passed_submissions{$p->{clone_name}};

    # look up its clone
    my $clone = CXGN::Genomic::Clone->retrieve_from_parsed_name($p);
    unless($clone) {
      warn "WARNING: no clone found for filename $p->{basename}";
      next;
    }

    # look up its chromosome assignment
    my $chr_in_registry = $clone->chromosome_num;
    unless( defined $chr_in_registry ) {
      warn "WARNING: $tarfile no longer has a chromosome assignment in the clone registry! (clone id $clone)";
      next;
    }

    # if they're the same, fine, go on to the next one
    next if $chr_in_registry eq $p->{chr};

    # otherwise, open the tarball, change its chromosome number, and republish it
    eval {

      my $s = CXGN::TomatoGenome::BACSubmission->open( $tarfile );

      $s->chromosome_number( $chr_in_registry );

      local $opt{f} = 1; #< fake like we have the -f option so it
                         #doesn't wait for a size changing check to
                         #process each of these
      process_file( $s->new_tarfile );

    }; if( $EVAL_ERROR ) {
      warn "failed to change $p->{basename} chromosome number from $p->{chr} to $chr_in_registry:\n$EVAL_ERROR";
      next;
    }

    print "successfully updated and republished $p->{basename}, changing chromosome from $p->{chr} to $chr_in_registry\n";
  }
}


# returns a list of validation errors that are NOT allowable, which
# can vary depending on whether the -G or -M options are passed to
# this script
sub _unallowable_errors(@) {
  my @errs = @_;
  my %allowable = ( G => { map $_ => 1,
			   E_GB_ACC,
			   E_GB_REC,
			   E_GB_SEQ, 
			   E_GB_PHASE_1,
			   E_GB_PHASE_3,
			 },
		    M => { E_MULT_SEQS, 1},
		    V => { E_VEC, 1},
		  );

  return if $opt{F};

  #filter out allowable errors, if any
  foreach my $allow_opt (keys %allowable) {
    if($opt{$allow_opt}) {
      @errs = grep {!$allowable{$allow_opt}{$_}} @errs;
    }
  }
  return @errs;
}


#given a submission object, set its version() using
#the sequences that are already published on the site.
#if there is a published sequence and this sequence is different,
#then make this one as the next version
sub set_submission_sequence_version {
  my $submission = shift;

  #look up the current published version for this bac
  my ($published_version,$published_seqs) =
    _get_published_seqs($submission->clone_object);

  if(DEBUG) {
    local $Data::Dumper::Maxdepth = 2;
    dprint "got published version $published_version, seqs ",Dumper($published_seqs);
  }

  unless( $published_version ) {
    $submission->version(1);
    return;
  }

  my @submitted_sequences = $submission->vector_screened_sequences;

  if( _same_sequences(\@submitted_sequences,$published_seqs) ) {
    print $submission->bac_name." sequence has not changed from version $published_version, not incrementing.\n";
    $submission->version($published_version);
  } else {
    my $new_version = $published_version + 1;
    print $submission->bac_name." sequence has changed since last version $published_version, new version is $new_version.\n";
    $submission->version($new_version);
  }
}


#take two arrays of Bio::Seq objects, return 1
#if they are the same sequences IN THE SAME ORDER
#otherwise, the fragment numbers would be messed up
#if they are the same sequences, but not ordered the
#same, they are considered a different set of sequences
#TODO: ask Lukas what he thinks about this
sub _same_sequences {
  my ($seqs1,$seqs2) = @_;

  return 0 unless @$seqs1 == @$seqs2;

  #well, now we know there are the same
  #number of seqs...

  foreach my $index (0..(@$seqs1-1)) {
    $seqs1->[$index]->seq eq $seqs2->[$index]->seq
      or return 0;
  }

  #if we haven't returned by now, they must be the same seqs
  return 1;
}

#get parsed idents and Bio::Seq objects for a given bac
sub _get_published_seqs {
  my $clone = shift;

  my $ident_index = _identifiers_index()
    or return ();

  my $published_seqs = $ident_index->{$clone->clone_id};
  dprint "111 for $clone, got from index:\n".Dumper $published_seqs;
  unless($published_seqs && @$published_seqs) {
    dprint "no clone $clone in identifiers index, index is:\n",Dumper $ident_index;
    return ();
  }

  #can take the first one, because we know our index will be consistent
  my @seqs =  map {
    my $ident = $_->{match};
    local $_; #the index thing apparently modifies $_.  protect us from that
    _seqs_index()->fetch( $ident )
      or die "$ident found in identifiers index, but not seqs index, parsed record is ".Dumper($ident);
  } @$published_seqs;
  dprint "222 for $clone, got from index:\n".Dumper $published_seqs;

  return ($published_seqs->[0]->{version},\@seqs);
}

#make an index by clone ID of the currently published sequence identifiers
#(in bacs.seq).  return this index as a hashref, as:
#{  clone_id => [ parsed identifier, parsed identifier, ...], ... }
memoize('_identifiers_index');
sub _identifiers_index {

  #find the name of the currently published version of bacs.seq
  my $allseqs_published_as = published_as(aggregate_filename('all_seqs'))
    or return();
  my $bacs_seqs_filename = $allseqs_published_as->{fullpath}
    or return {}; #return an empty index if there's no published aggregates file

  #this is the index we're building
  my %identifiers_index;

  #parse the identifiers in the published sequence file to find their
  #versions, all the fragments, etc., and build the index of them
  open my $seqfile, $bacs_seqs_filename
    or die "Could not open published bac sequences '$bacs_seqs_filename' for reading: $!";
  while(<$seqfile>) {
    if(my ($identifier) = /^\s*>\s*(\S+)/) {
      my $parsed_ident = parse_clone_ident($identifier,'agi_bac_with_chrom','versioned_bac_seq')
	or die "Unparsable identifier ($identifier) found in published sequences file";
#      warn "got parsed ident ".Dumper $parsed_ident;
      my $clone = CXGN::Genomic::Clone->retrieve_from_parsed_name($parsed_ident)
	or die "No clone with identifier '$identifier' found in database";

      $identifiers_index{$clone->clone_id} ||= [];
      push @{$identifiers_index{$clone->clone_id}},$parsed_ident;
    }
  }

#  warn 'returning idents index ',Dumper \%identifiers_index;

  #return ref to the index
  return \%identifiers_index;
}

#read the bac sequences and index them by sequence identifier
#Memoized, so this is only done once per run of this script.
#Subsequent runs return the same index object
memoize('_seqs_index');
sub _seqs_index {

  #find the name of the currently published version of bacs.seq
  my $bacs_seqs_filename = published_as(aggregate_filename('all_seqs'))->{fullpath}
    or return {}; #return an empty index if there's no published aggregates file

  #make an index of that sequence file
  my (undef,$tempfile) = tempfile(UNLINK=>1);
  my $i = Bio::Index::Fasta ->new(-filename => $tempfile, -write_flag => 1)
    or die "Could not make bac seqs index from '$bacs_seqs_filename' into temp file $tempfile (error 1): $!";
  $i->make_index( $bacs_seqs_filename )
    or die "Could not make bac seqs index from '$bacs_seqs_filename' into temp file $tempfile (error 2): $!";

  #and return the index
  return $i;
}


