#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use Bio::SeqIO;
use Bio::Seq::Quality;

use CXGN::TomatoGenome::Config;
use CXGN::TomatoGenome::BACSubmission qw/E_MULT_SEQS E_GB_ACC E_GB_REC E_GB_SEQ/;
use CXGN::Tools::List qw/str_in list_join flatten/;
use CXGN::Tools::Run;

our $vhost = CXGN::TomatoGenome::Config->new;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] orig_file new_filename

  Given an unfinished BAC in the old multi-sequence format, reformat
  it as a single sequence stuck together with N's.

  Options:

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();

my ($subfile,$targetfile) = @ARGV;
$subfile && -f $subfile or usage;
-f $targetfile and die "Cowardly refusing to overwrite $targetfile\n";

my $sub = CXGN::TomatoGenome::BACSubmission->open($subfile);

#check for validation errors that we can't ignore
my @errors = $sub->validation_errors;
my @allowed_errors = ( E_MULT_SEQS, E_GB_ACC, E_GB_REC, E_GB_SEQ );
grep !str_in($_,@allowed_errors),@errors
  and die "Validation failed in a way that is not tolerated by this script:\n".$sub->validation_text;
#check that we have the mult seqs error
str_in(E_MULT_SEQS,@errors)
  or die "Submission $subfile does not contain multiple seqs.\n";

#now get all our sequences and tack them together
my $seqfile = $sub->sequences_file;
my $bigseq = do {
  my $seq_in = Bio::SeqIO->new( -file => $seqfile, -format => 'fasta');
  my @seqs;
  while (my $s = $seq_in->next_seq) {
    push @seqs, $s->seq;
  }
  Bio::PrimarySeq->new( -display_id => $sub->sequence_identifier, -seq => join('N'x100,@seqs));
};
#now write out the concatenated seq
Bio::SeqIO->new( -file => ">$seqfile", -format => 'fasta')->write_seq($bigseq);

my $qualfile = $sub->qual_file;
if (-f $qualfile) {
  my $bigqual = do {
    my $qual_in = Bio::SeqIO->new( -file => $qualfile, -format => 'qual' );
    my @quals;
    while ( my $q = $qual_in->next_seq ) {
      push @quals, $q->qual;
    }
    Bio::Seq::Quality->new( -id => $sub->sequence_identifier, -qual => [flatten( list_join [(0)x100], @quals )]);
  };
  Bio::SeqIO->new( -file => ">$qualfile", -format => 'qual')->write_seq($bigqual);
}

#make sure that it no longer has multi-sequence errors
die "Sequence reformat failed, still has multiple sequence validation error:\n".$sub->validation_text
  if str_in(E_MULT_SEQS,$sub->validation_errors);

CXGN::Tools::Run->run( 'cp', $sub->new_tarfile, $targetfile );


