#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Basename;
use File::Spec;
use File::Copy;

#use Data::Dumper;

use Bio::SeqIO;

use CXGN::TomatoGenome::BACSubmission;
use CXGN::Tools::Run;


sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script -o <file_or_dir> submission_file submission_file ...

  Create a properly-formatted GenBank submission file for a sequenced
  BAC using the NCBI fa2htgs program, outputting it to the given file
  or directory.  If a directory is given as the -o option, then files
  will created in that directory named for each BAC.

  Options:

  -o <file or dir>
    output the formatted submission to this file, or for multiple
    submissions, this directory

  -c <center_name>
    force the name of the sequencing center to use for generating the
    submission.  By default, the sequencing center is figured out
    based on the chromosome of the submitted BAC.

  -t <template_file>
    force the use of the given Sequin-generated submission template
    file.  Normally, the Sequin template in
    $FindBin::RealBin/templates/<seq_center>.sqn is used.

  -p <phase>
    force an HTGS phase on each of the submissions.  By default,
    guesses for each tarball based on whether the sequence looks
    finished.


  PLEASE NOTE:

    If some BACs have already been submitted to GenBank, make sure
    their submission tarballs contain the correct gbacc.txt file
    before running this script on them.  This is required for the
    submission to be properly flagged as an update.

EOU
}

our %opt;
getopts('o:c:t:p:',\%opt) or usage();
$opt{o} or usage;
@ARGV or usage;

foreach my $bsub_file (@ARGV) {
  #figure out the output file to hold the genbank submission

  eval {
    my $sub = CXGN::TomatoGenome::BACSubmission->open($bsub_file);
    my $gbfile = $sub->genbank_submission_file($opt{p},$opt{c},$opt{t});
    my $bn = $sub->clone_object->clone_name_with_chromosome;
    my $outfile = -d $opt{o} ? File::Spec->catfile($opt{o},"$bn.htgs") : $opt{o};
    copy($gbfile,$outfile);
  }; if($EVAL_ERROR) {
    warn "Failed to make a GenBank file for $bsub_file:\n  $EVAL_ERROR\n";
  }
}
