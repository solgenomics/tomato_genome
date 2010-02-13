#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use Hash::Util qw/lock_hash/;

#use Data::Dumper;
use File::Basename;
use File::Spec;
use List::Util qw/sum/;

use Bio::FeatureIO;
use Bio::SeqIO;

use CXGN::DB::Connection;

use CXGN::IndexedLog;

use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident assemble_clone_ident/;
use CXGN::Publish qw/move_or_print copy_or_print/;

use CXGN::TomatoGenome::BACSubmission ':errors';
use CXGN::TomatoGenome::BACPublish qw/publisher parse_filename/;
use CXGN::TomatoGenome::Config;

use CXGN::Tools::File qw/file_contents/;

##### DEFAULTS #######

my $cfg = CXGN::TomatoGenome::Config->load_locked;
our $country_topdir = $cfg->{'country_uploads_path'};
our $bac_pipeline_dir = $cfg->{'bac_pipeline_dir'};

######################

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options]  country_name country_name ...

  Script to run preprocessing steps on a BAC submission before it is
  published.  Meant to be run by cron on several country directories.

  Finds all submission tarballs in the country_dir/pre-upload
  directories and preprocesses them to prepare them for submitting to
  GenBank and publishing on SGN.  Processed submissions are moved into
  pre-upload-processing while being processed, and then into
  pre-upload-processed when they've been successfully processed.
  Also, successfully processed files are _copied_ into upload/ for
  processing by the main submission pipeline.

  Steps in the preprocessing:

    - vector trim: trim the ends of the input sequence to remove
                   vector sequence

    - genbank submission: generate a .htgs file for submitting to GenBank
                          and put it in the <country>/genbank_submissions/new
                          directory

  Options:

  -V   skip the vector trimming step

  -G   skip the GenBank submission step

  -v   be verbose

  -d <dir>
    set top-level country upload directory
    Default: $country_topdir

  -x  dry run.  don't actually change anything (implies -v)

  -r  attempt to repair BAC submissions if they don't validate

  -c <center>
      force a sequencing center name to be used for all the BACs

  -t <file>
      force a Sequin template to use for all the BACs

  -P  vector trim even if the detected vector regions don't contain T7
      and SP6 promoter seqs

EOU
}

our %opt;
check_executables();
getopts('VGvxrc:d:t:P',\%opt) or usage();
$CXGN::Publish::dry_run = 1 if $opt{x};
$CXGN::Publish::print_ops = 1 if $opt{x} || $opt{v};
$country_topdir = $opt{d} if defined $opt{d};
our $local_genbank_top_dir = File::Spec->catdir( $bac_pipeline_dir, $cfg->{'bac_genbank_dir'} );
-d $local_genbank_top_dir or mkdir $local_genbank_top_dir; #< try to make it if it's not there

sub vprint(@) {
  print @_ if $opt{v} || $opt{x};
}

my @country_dirs = map {File::Spec->catdir($country_topdir,$_)} @ARGV
  or usage;

foreach my $dir (@country_dirs) {
  -d $dir or die "No such directory '$dir'\n";
}

my %allowable_errors = map {$_=>1} ( E_GB_ACC,  #< no genbank accession
				     E_GB_REC,  #< malformed genbank record
				     E_GB_SEQ,  #< genbank seq doesn't match this seq
				   );


########
### generate genbank files for all submissions that need preprocessing
########
foreach my $country_dir (@country_dirs) {
  my @subfiles = glob(File::Spec->catfile($country_dir,'pre-upload','*.tar.gz'));
  #check for duplicate submissions for BACs
  { my %seen;
    foreach my $sf (@subfiles) {
      my $p = parse_filename($sf);
      my $cn = $p->{clone_name};
      if($seen{$cn}) {
	die "duplicate submissions found for $cn: $seen{$cn} and $sf, aborting\n";
      }
      $seen{$cn} = $sf;
    }
  }
  foreach my $subfile (@subfiles) {
    my ($bn,$dirname) = fileparse($subfile);

    my $upload_loc = $subfile;
    my ($done_loc) = map {
      my $tag = $_;
      my $d = $dirname;
      $d =~ s!/$!-$tag!;
      File::Spec->catfile($d,$bn);
    } qw/processed/;

    eval {
      my $sub = CXGN::TomatoGenome::BACSubmission->open($upload_loc);

      #if -U is given, stick the unfinished sequences together with a certain number of Ns
      if( $sub->sequences_count > 1 ) {
	vprint "concatenating sequences for multi-sequence submission (-U option given)\n";

	my $seqs = Bio::SeqIO->new(-file => $sub->sequences_file, -format => 'fasta');
	my $firstseq = $seqs->next_seq;

	#fix the identifier on the first sequence
	my $f_ident = parse_clone_ident($firstseq->primary_id,'versioned_bac_seq');
	$f_ident  or die "could not parse '".$firstseq->primary_id."'\n";
	delete $f_ident->{fragment};
	$firstseq->primary_id(assemble_clone_ident(versioned_bac_seq => $f_ident) or die 'sanity check failed');

	#now stick all the others to it
	while( my $s = $seqs->next_seq ) {
	  $firstseq->seq( $firstseq->seq.('N' x 100).$s->seq );
	}
	$seqs = undef;

	#write out a new sequences file with a single seq
	Bio::SeqIO->new(-file => '>'.$sub->sequences_file,
			-format => 'fasta',
		       )->write_seq($firstseq);


	#if there is a qual file, fix that too
	if( -f $sub->qual_file ) {
	  die "FIX ME, qual file processing for fragmented sequences not yet implemented";
	}
      }

      #attempt to repair, then skip this one if it has an error that
      #is not one of the allowable ones
      $sub->repair if($opt{r} && grep {!$allowable_errors{$_}} $sub->validation_errors);
      if(grep {!$allowable_errors{$_}} $sub->validation_errors) {
	die $sub->validation_text;
      }

      #sub_changed keeps track of whether we've actually done anything to this submission
      my $sub_changed = preprocess($sub);

      #look up any existing genbank accession for this submission,
      #setting it inside this submission if necessary
      if(!$sub->genbank_accession && $sub->clone_object->genbank_accession) {
	$sub_changed = 1;
	$sub->genbank_accession($sub->clone_object->genbank_accession);
      }

      #make a genbank submission file for this
      my $outgoing_dir = File::Spec->catdir($local_genbank_top_dir,'outgoing');
      -d $outgoing_dir or mkdir $outgoing_dir; #< try to make it if not there
      publisher()->publish(['cp',$sub->genbank_submission_file($opt{c},$opt{t}),File::Spec->catfile($outgoing_dir,$sub->clone_object->clone_name_with_chromosome.'.htgs')]);

      #now that we're finished with it:
      # - move it into the pre-upload-processed subdir of this country dir
      # - copy it into the upload/ subdir of this country dir, where it will sit and be stuck in the pipeline until we add an accession to it, and genbank publishes the record
      move_or_print($upload_loc,$done_loc)
	or die "could not move $upload_loc => $done_loc: $!";
      my $upload_dir = File::Spec->catdir($country_dir,'upload');
      -d $upload_dir or mkdir $upload_dir; #< try to make it if not there
      my $next_stage_target = File::Spec->catfile($upload_dir,$sub->clone_object->clone_name_with_chromosome.'.tar.gz');
      if( $sub_changed ) {
	move_or_print($sub->new_tarfile, $next_stage_target)
	  or die "could not move new tarfile to $next_stage_target: $!";
      } else {
	copy_or_print($done_loc, $next_stage_target)
	  or die "could not copy $done_loc to $next_stage_target: $!";
      }
    }; if($EVAL_ERROR) {
      warn $EVAL_ERROR;
    }
  }
}

########
### upload to genbank any .htgs files that have not yet been uploaded
########

my %ncbi = (  server     => 'ftp-private.ncbi.nlm.nih.gov',
	      user       => 'cornell',
	      password   => '4dluvO$$',
	      upload_dir => 'SEQSUBMIT',
	      report_dir => 'REPORT',
	   );
lock_hash %ncbi;

my $dbh = CXGN::DB::Connection->new;
my $log = CXGN::IndexedLog->open( DB => $dbh,
				  $cfg->{'genbank_upload_log'}
				);
foreach my $country_dir (@country_dirs) {
  my @htgs_files = glob(File::Spec->catfile($local_genbank_top_dir,'outgoing','*.htgs'));

  foreach my $htgs_file (@htgs_files) {
    if( $log->lookup(content => "UPLOADED_TO_GENBANK $htgs_file") ) {
      vprint("not uploading $htgs_file to genbank, already logged as done.\n");
      next;
    }

    #upload it to ncbi
    eval {
      unless( $opt{x} ) {
	vprint("uploading $htgs_file to ncbi...\n");
	my $put = CXGN::Tools::Run->run('ncftpput',
					$ncbi{user} ? ( -u => $ncbi{user} ) : (),
					$ncbi{password} ? (-p => $ncbi{password} ) : (),
					$ncbi{server},
					$ncbi{upload_dir},
					$htgs_file,
				       );
      } else {
	vprint("dry run, skipping upload of $htgs_file to genbank\n");
      }
    }; if($EVAL_ERROR) {
      warn "$htgs_file upload failed: $EVAL_ERROR\n";
      $log->append("Genbank upload failed for $htgs_file: $EVAL_ERROR") unless $opt{x};
    } else {
      $log->append("UPLOADED_TO_GENBANK $htgs_file") unless $opt{x};
    }
  }
}

######
# get any new reports from ncbi.  'wget -nc' will not re-download files
# that are already present locallly.
######

my $local_report_dir = File::Spec->catfile( $local_genbank_top_dir, 'reports_from_ncbi' );
-d $local_report_dir or mkdir $local_report_dir; #< make it if possible
unless( $opt{x} ) {
  my $dl = CXGN::Tools::Run->run('wget',
				 '-nc',
				 $ncbi{user}     ? ( "--ftp-user=$ncbi{user}" ) : (),
				 $ncbi{password} ? ( "--ftp-password=$ncbi{password}" ) : (),
				 "ftp://$ncbi{server}/$ncbi{report_dir}/*.ac4htgs",
				 {
				  working_dir => $local_report_dir,
				 },
				);
} else {
  vprint("dry run, skipping download sync from genbank\n");
}


#######
# now add gb accessions to any files in country upload/ dirs that match report files
#######

# look for submissions that are (still) in the beginning of the main
# pipeline
my @uploaded_submissions = glob( File::Spec->catfile( $country_topdir, '*', 'upload', 'C*.tar.gz') );

foreach my $subfile (@uploaded_submissions) {
  my $p = parse_filename($subfile); #< from CXGN::TomatoGenome::BACPublish

  vprint("looking for report file for $p->{basename}...\n");

  # do we have any ac4htgs files for this uploaded submission?
  my @reportfiles = glob( File::Spec->catfile( $local_report_dir, "*.$p->{clone_name}.*.ac4htgs" ) )
    or next; #< if not, skip to the next submission

  # sort the report files by their ncbi upload date
  @reportfiles = sort {
    # because of how ncbi formats these report file names, the first
    # set of numbers in each filename will be the date formatted like
    # YYYYMMDD
    my ($a_date) = basename($a) =~ /(\d+)/;
    my ($b_date) = basename($b) =~ /(\d+)/;

    $a_date <=> $b_date
  } @reportfiles;


  # take the most recent report file from the date-sorted array
  my $most_recent_reportfile = pop @reportfiles;

  vprint("report file $most_recent_reportfile found, parsing.\n");

  my $accession = get_unversioned_accession_from_ac4htgs_file( $most_recent_reportfile );

  vprint("parsed accession '$accession'.\n");
  my $sub = CXGN::TomatoGenome::BACSubmission->open($subfile);

  # save the state of the submission before we start messing with it
  my $old_acc = $sub->genbank_accession || '';
  my $old_err_cnt  = $sub->validation_errors;
  my $old_err_txt = $sub->validation_text;

  if( $old_acc =~ /^$accession/ ) {
    vprint("old accession is the same as the new one, skipping.\n");
    next;
  }

  vprint("setting new accession '$accession' in $subfile (old accession was '$old_acc')...\n");

  # now set the accession
  $sub->genbank_accession($accession);

  # did we reduce its validation errors by setting this accession? we
  # don't want to make this modification if we made things worse
  my $new_err_cnt = $sub->validation_errors;
  my $new_err_txt = $sub->validation_text;

  if( $new_err_txt ) {
    warn "warning: still had errors after setting accession $accession on $p->{basename}:\n";
    warn $new_err_txt;
  }

  if ( $old_err_cnt < $new_err_cnt ) {
    $old_err_txt ||= "none\n";
    $old_acc ||= 'none';
    warn( <<EOF );
old errors in $p->{basename}, with containing old accession ($old_acc):
$old_err_txt
new errors are more numerous than old errors.  maybe '$accession' is not the right accession to be setting?\n";
EOF
    next;
  }

  # now finally replace the uploaded submission with a new one
  move_or_print( $sub->new_tarfile, $subfile )
    or die "$! replacing submission file $subfile";
}


########################## SUBROUTINES #############################

#run preprocessing steps on this submission, currently just vector screening
sub preprocess {
  my ($sub) = @_;
  if(vector_trim($sub)
     # OR-in other altering procedures here
    ) {
    return $sub->new_tarfile;
  }
  return;
}

#trim the sequence inside this submission if it contains two big hunks of vector.
#returns 1 if it actually did any trimming
sub vector_trim {
  my ($sub) = @_;
  #skip unless there's some vector sequence in it
  return unless file_contents($sub->vector_screened_sequences_file) =~ /X+/i;

  my (undef,$vector_gff3) = $sub->analyze_with('Cross_match_vector');

  my $seq = Bio::SeqIO->new( -file       => $sub->sequences_file,
			     -format     => 'fasta',
			   )->next_seq
			     or confess 'no seq!';

  #attach the annotations to this seq
  my $vec_feats = Bio::FeatureIO->new( -file    => $vector_gff3,
				       -format  => 'gff',
				       -version => 3,
				     );
  my @vec_feats;
  my $vec_feat_count = 0;
  while(my $f = $vec_feats->next_feature) {
    $vec_feat_count++;
#    warn "got vec feat (".$f->start.",".$f->end.")\n";
    $f->attach_seq($seq);
    push @vec_feats,$f;
  }
  return unless $vec_feat_count;

  if( $vec_feat_count >= 2 ) {
     #try to find the SP6 and T7 promoter sites that bracket the insert sequence.
     my $t7_seq  = 'TAATACGACTCACTATAGGG';
     my $sp6_seq = 'ATTTAGGTGACACTATAG';
     sub revcomp($) { my $seq = reverse $_[0]; $seq =~ tr/ATCG/TAGC/; $seq}
     sub contains { index($_[0],$_[1]) != -1 }

#     print revcomp $_,"\n" foreach $t7_seq,$sp6_seq;

     #find features that are large, and lie no more than <ends_window> from the
     #end of a sequence
     my $ends_window = 1_000;
     my $max_trimmed_bases = 10_000; #< maximum number of trimmed-off bases we consider OK

     my $first_bit   = Bio::Range->new(-start => 1,                         -end => $ends_window);
     my $last_bit    = Bio::Range->new(-start => $seq->length-$ends_window, -end => $seq->length);
     my $whole_thing = Bio::Range->new(-start => 1,                         -end => $seq->length);

     my $left_trim;
     my $right_trim;
     foreach my $f (@vec_feats) {
       next unless $f->length >= 2000;
       if($f->overlaps($first_bit)) {
	 $left_trim = $left_trim ? $f->union($left_trim) : $f;
       } elsif($f->overlaps($last_bit)) {
	 $right_trim = $right_trim ? $f->union($right_trim) : $f;
       }
     }

     #return unless we have viable candidates for trimming
     return unless $left_trim && $right_trim;

     $left_trim->start(1);
     $right_trim->end($seq->length);

     #if it's too much vector, die
     my $trimmed_bases = $left_trim->length + $right_trim->length;
     $trimmed_bases <= $max_trimmed_bases
       or die "This sequence seems to have $trimmed_bases that need trimming off.  That's more than my maximum threshold of $max_trimmed_bases";

     sub r2str { '('.$_[0]->start.','.$_[0]->end.')'}
     vprint r2str($left_trim)." and ".r2str($right_trim)." from sequence (1,".$seq->length.") were masked with vector.  Do they contain T7 and SP6 promoters?\n";

     #the candidate vector regions must contain the T7 and SP6 promoter sequences
     $left_trim->attach_seq($seq);
     $right_trim->attach_seq($seq);
     unless( $opt{P}
	     or contains($left_trim->seq->seq,$t7_seq) && contains($right_trim->seq->seq,revcomp($sp6_seq))
	     or contains($left_trim->seq->seq,$sp6_seq) && contains($right_trim->seq->seq,revcomp($t7_seq))
	   ) {
       vprint "No promoter sequences found, not trimming\n";
       return;
     }

     vprint "Promoter sequences found, trimming sequence.\n";

     #now trim the sequence
     $seq->seq( $seq->subseq( $left_trim->end+1, $right_trim->start-1 ) );

     #write out the new sequence file
     Bio::SeqIO->new( -file => '>'.$sub->sequences_file,
		      -format => 'fasta',
		    )->write_seq($seq);

     return 1;
   }
   return;
}


# given a *.ac4htgs file, open it and parse the accession out of it,
# remove the version from the accession, and return a string
# containing the unversioned accession
sub get_unversioned_accession_from_ac4htgs_file {
  my ($reportfile) = @_;

  #now parse the accession out of the report file
  open my $rfh, $reportfile or die "$! opening ncbi report file '$reportfile'";
  my $accession;
  while(<$rfh>) {
    if(/accession\s*:/i) {
      ($accession) = /accession\s*:\s+([\w\.]+)/i;
      last;
    }
  }

  #make sure we got a sensible-looking accession
  $accession =~ /^[A-Z]{2}\d+\.\d+$/
    or die "parsed out invalid accession '$accession' from report file '$reportfile'";


  #take the version number off of the accession
  $accession =~ s/\.d+$//;

  return $accession;
}


# check that we have all the executables that we need
sub check_executables {
  my @required = qw/ncftpput wget/;
  my $err = 0;
  foreach (@required) {
    unless(`which $_`) {
      $err = 1;
      warn "this script requires ncftpput executable in your path\n";
    }
  }

  $err and die "not all required executables present, aborting.\n";
}
