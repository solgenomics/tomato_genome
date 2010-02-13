#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use CXGN::TomatoGenome::BACSubmission;

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script new_chromosome_number submission_tarball
  Example: $FindBin::Script 5 C04HBa0001A01.tar.gz

  Script to change the chromosome number in a BAC submisssion tarball
  and then run bsub_publish_bac_seq.pl on the tarball with changed
  chromosome number.

  Options:

    -d dir
       publishing directory to republish to, passed along to
       bsub_publish_bac_seq.pl

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('d:',\%opt) or usage();


my ($new_chrom,$tarfile) = @ARGV;

my $s = CXGN::TomatoGenome::BACSubmission->open($tarfile);

$s->chromosome_number($new_chrom);

CXGN::Tools::Run->run( 'bsub_publish_bac_seq.pl',
		       '-f',
		       $opt{d} ? (-d => $opt{d}) : (),
		       $s->new_tarfile
		      );

