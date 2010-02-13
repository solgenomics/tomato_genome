#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Spec;

use CXGN::TomatoGenome::Config;

use CXGN::TomatoGenome::BACPublish qw/ find_submissions /;

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script

  Search upload dirs for BAC submissions and process them.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();

my $cfg = CXGN::TomatoGenome::Config->load_locked;
my $ftpdir = File::Spec->catdir( @{$cfg}{ 'ftpsite_root',
                                          'bac_publish_subdir'
                                        },
                               );

# look for submissions
my @submission_files = find_submissions($opt{d});

exit unless @submission_files;

print map "$_\n",
  'cron job found new submissions:',
  map "   $_", @submission_files;

print "running bsub_publish_bac_seq.pl on them:\n";

exec 'bsub_publish_bac_seq.pl', '-c', '-S', -d => $ftpdir, @submission_files;
