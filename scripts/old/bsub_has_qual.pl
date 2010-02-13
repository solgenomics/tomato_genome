#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use CXGN::TomatoGenome::BACSubmission;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] <submission> <submission> ...

  Tiny script, prints "<seqname> has qual" if the given submission(s)
  have main qual files included.

  Options:

    -n   if passed, print '<seqname> does not have qual' if the
         submission doesn't have a qual file.  off by default.

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('n',\%opt) or usage();

foreach my $file (@ARGV) {
  my $sub = CXGN::TomatoGenome::BACSubmission->open_stripped($file)
    or die "$! opening $file";

  if( -f $sub->qual_file ) {
    print $sub->sequence_identifier." has qual\n";
  }
  elsif( $opt{n} ) {
    print $sub->sequence_identifier." does not have qual\n";
  }
}

