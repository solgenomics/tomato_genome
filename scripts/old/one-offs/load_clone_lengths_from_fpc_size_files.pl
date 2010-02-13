#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident/;
use CXGN::Genomic::Clone;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script size_file size_file ...

  Takes a list of FPC band size files (usually called
  <clone_name>.sizes) on the command line, parses each one, adds up
  the band sizes for each clone to get the clone's estimated length,
  then loads the estimated length in the genomic.clone table.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();


my @sizefiles = @ARGV;
foreach (@sizefiles) {
  -r or die "can't open $_ for reading";
}


#make sure all the clone names in each size file are parsable
foreach my $sizefile (@sizefiles) {
  open my $sfh, $sizefile or die "$! reading $_";
  while(<$sfh>) {
    next unless /[^-\d\s]/; #skip lines that don't contain clone names
    my ($clone_name) = split;
    my $p = parse_clone_ident($clone_name)
      or die "cannot parse clone name identifier $clone_name.  either it's corrupted, or you just need to implement support for it in CXGN::Genomic::CloneIdentifiers\n";
  }
  close $sfh;
}

#now load all the lengths
foreach my $sizefile (@sizefiles) {
  open my $sfh, $sizefile or die "$! reading $_";
  while(<$sfh>) {
    chomp;
    if( /[^-\d\s]/ ) { #clone name line
      my ($clone_name) = split;
      my $clone = CXGN::Genomic::Clone->retrieve_from_clone_name($clone_name)
	or die "cannot look up clone from name $clone_name.  either it's corrupted, or you just need to implement support for it in CXGN::Genomic::CloneIdentifiers\n";
      my $length = 0;
      while(<$sfh>) {
	chomp;
	last if $_ == -1;
	$length += $_;
      }
      if( $clone->estimated_length ) {
	warn "WARNING: $clone_name already has estimated length ".$clone->estimated_length."\n";
	next;
      } else {
	#warn "setting $clone_name length $length\n";
	$clone->estimated_length($length);
	$clone->update;
      }
    } else {
      die 'parse error';
    }
  }
  close $sfh;
}


CXGN::Genomic::Clone->db_Main->commit;
