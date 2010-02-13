#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use CXGN::Tools::List qw/str_in/;
use CXGN::Tools::Identifiers qw/identifier_namespace/;

use CXGN::BioTools::AGP qw/agp_to_seq agp_parse agp_contig_seq agp_contigs /;

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script agp_file

  Given an AGP file on the command line, assembles either contig
  sequences (the default), or a large N-gapped pseudomolecule sequence for it.
  Currently knows how to fetch sequences for the following
  CXGN::Tools::Identifiers namespaces:
     bac_sequence

  Options:

  -P make one big pseudomolecule sequence instead of contig sequences

  -g <number>
     set a fixed length override for gaps.  Defaults to 0, meaning
     make gaps as large as called for in the AGP

  -n <letter>
     if -P, fill gaps with given character
     default N

  -l
     if set, output lower-case sequence
     default off

  -d <defline string>
     set the fasta definition line to the given string

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('g:n:lP',\%opt) or usage();
$opt{g} ||= 0;
$opt{g} >= 0 or usage("-g must be >= 0");
$opt{n} ||= 'N';
length($opt{n}) == 1 or usage("-n must be a single letter");
$opt{d} && index($opt{d},"\n") != -1 and usage("-d defline cannot contain a newline");

usage('invalid arguments') unless @ARGV == 1;

-r $ARGV[0] or die "cannot open $ARGV[0] for reading\n";

if( $opt{P} ) {
  my ($ident,$seq) = agp_to_seq($ARGV[0],
				fetch_bac_sequence => \&fetch_bac_sequence,
				gap_length  => $opt{g},
				no_seq_char => $opt{n},
				lowercase   => $opt{l},
			       );

  print '>'.$ident;
  print ' '.$opt{d} if $opt{d};
  print "\n";
  print $seq,"\n";

} else {
  my @contigs = agp_contigs(agp_parse($ARGV[0]));

  my $ctr;
  foreach my $ctg (@contigs) {
    my $seq = agp_contig_seq( $ctg, fetch_bac_sequence => \&fetch_bac_sequence, lowercase => $opt{l} );
    print '>'.$ARGV[0].'.contig'.++$ctr;
    print ' '.$opt{d} if $opt{d};
    print "\n";
    print $seq,"\n";
  }
}

#### SEQUENCE FETCHING ROUTINES ####
# each function is named type_fetch
# where type is the string returned by a call to CXGN::Tools::Identifiers::identifier_namespace

sub fetch_bac_sequence {
  my ($ident) = @_;
  return _get_chado_seq($ident);
}

#fetch a sequence by name from chado
sub _get_chado_seq {
  my ($ident) = @_;
  our $dbh ||= CXGN::DB::Connection->new;
  our $chado_feature_seq_q ||= $dbh->prepare('select residues from feature where name = ?');
  $chado_feature_seq_q->execute($ident)
    or return;
  my ($seq) = $chado_feature_seq_q->fetchrow_array;
  $chado_feature_seq_q->finish;
  return $seq;
}
