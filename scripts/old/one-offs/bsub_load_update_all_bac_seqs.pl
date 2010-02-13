#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use Bio::SeqIO;

use CXGN::TomatoGenome::BACPublish qw/aggregate_filename publisher/;

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script -v <publishing directory>

  Update the residues column of the chado feature table using the
  sequences from the 'all repeatmasked seqs' file in the given bac
  publishing directory.

  Options:

    -v  be verbose

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('v',\%opt) or usage();

sub vprint(@) {
  print @_ if $opt{v};
}
-d $ARGV[0] or usage;
@ARGV == 1 or usage;

my $seqfile = aggregate_filename('all_seqs',$ARGV[0]);
$seqfile = publisher->published_as($seqfile)
  or die "no repeatmasked bac seqs file found, is the given publishing dir correct?";
$seqfile = $seqfile->{fullpath};

my $seq_in = Bio::SeqIO->new(-file => $seqfile, -format => 'fasta'); #< dies on error

my $dbh = CXGN::DB::Connection->new({ dbargs => {AutoCommit=>0}});
my $check_q  = $dbh->prepare('select count(*) from feature where name=?');
my $update_q = $dbh->prepare('update feature set residues=?,uniquename=name,seqlen=? where name=?');
while(my $seq = $seq_in->next_seq) {
  my $name = $seq->display_id;
  $check_q->execute($name);
  my ($cnt) = $check_q->fetchrow_array;
  unless($cnt == 1) {
    warn "$cnt features with name $name, skipping\n";
    next;
  }

  my $seqlen = $seq->length;
  vprint "updating $name with $seqlen bases\n";
  $update_q->execute($seq->seq,$seqlen,$name);
}

$dbh->commit;
