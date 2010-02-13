#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use CXGN::TomatoGenome::BACPublish qw/parse_filename publishing_locations publisher/;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [ pub_dir ]

  Look in the given BAC publishing directory, find obsolete files, and
  unpublish them.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();


my ($pub_dir) = @ARGV;

my @files = grep -f, glob( "$pub_dir/chr*/*finished/C{0,1}*.tar.gz" );

my $publisher = publisher();

foreach my $file (@files) {
  next unless -f $file; #< might have already been deleted
  print "checking obsoletes for '$file'...\n";
  my $p = parse_filename( $file )
    or die "could not parse $file";

  my $pub = publishing_locations( $pub_dir, $p->{seq_name}, $p->{finished} );

  foreach my $obs (map @$_, values %{$pub->{obsolete}}) {
    if( my $published = $publisher->published_as($obs) ) {
      print "removing obsolete file '$published->{fullpath}'\n";
      $publisher->publish([rm => $published->{fullpath}]);
    }
  }

}
