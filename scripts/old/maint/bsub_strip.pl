#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Copy;

#use Data::Dumper;
use CXGN::TomatoGenome::BACSubmission;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script -o outfile  submission_file

  Script to strip all non-critical data out of a BAC submission
  tarball, making a new, much smaller one, for the purposes of testing
  BAC submission handling scripts.

  Options:

  -o outfile
    File to copy the new tarball to.

EOU
}

our %opt;
getopts('o:',\%opt) or usage();

my $sub = CXGN::TomatoGenome::BACSubmission->open(shift @ARGV);

my $subdir = $sub->main_submission_dir;
-d $subdir or die "No main submission directory '$subdir' found inside tarball";

#delete all subdirectories
system "rm -rf $subdir/*/";
die "rm -rf $subdir/*/ failed: $!" if $CHILD_ERROR;

my $tarball = $sub->new_tarfile;
copy($tarball,$opt{o})
  or die "could not copy $tarball -> $opt{o}: $!\n";
