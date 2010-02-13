#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use Memoize;

#use Data::Dumper;
use CXGN::Genomic::Clone;
use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident assemble_clone_ident/;
use CXGN::TomatoGenome::BACPublish qw/agp_file/;
use CXGN::BioTools::AGP qw/agp_parse/;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script bac_list

  Given a list of BACs as a file on the command line, or on STDIN, looks up the position
  of each one in the AGP files and prints it in CSV format on STDOUT.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();


print '"BAC name", "Chromosome", "Start Base","Start % from N", "End Base","End % from N"'."\n";

while( my $bacname = <> ) {
  chomp $bacname;

  my $c = CXGN::Genomic::Clone->retrieve_from_clone_name( $bacname )
    or next;

  my $chr = $c->chromosome_num
    or next;

  my $a = indexed_agp($chr);
  my $l = agp_line($c);

  my ($start,$sp,$end,$ep) = do {
    if( $l ) {
      my $s = sprintf('%0.2f',$l->{ostart}/$a->{max_base}*100);
      my $e = sprintf('%0.2f',$l->{oend}/$a->{max_base}*100);

      $l->{ostart},$s,$l->{oend},$e;
    } else {
      ('')x4
    }
  };

  print join ',',$c->clone_name, $chr, $start,$sp,$end,$ep;
  print "\n";

}

sub agp_line {
  my $c = shift;

  my $chr = $c->chromosome_num;
  my $name_with_chr = $c->clone_name_with_chromosome;

  my $a = indexed_agp($chr);
  return $a->{$name_with_chr};
}

sub indexed_agp {
  my ($chr) = @_;

  my $agp = agp_file($chr);

  my $p = agp_parse($agp);

  my %index;
  foreach my $line (@$p) {
    next unless $line->{ident};
    my $p = parse_clone_ident($line->{ident},'versioned_bac_seq')
      or die "could not parse '$line->{ident}'";
    my $n = assemble_clone_ident( agi_bac => $p );
    warn "WARNING: duplicate lines found for $n\n" if $index{$n};
    $index{ $n } = $line;
  }

  $index{max_base} = $p->[-1]->{oend};

  return \%index;
}
BEGIN { memoize('indexed_agp'); }
