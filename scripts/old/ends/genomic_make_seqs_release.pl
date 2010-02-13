#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Basename;
use File::Temp;

use CXGN::Genomic;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script  dirname

  This script make a BAC end sequence release, consisting of a bunch
  of sequence files, qual files, and formatted BLAST db files, in the
  given directory.

  Options:

   -L libname,libname,...

      list of case sensitive library shortnames to include.  Use
      'combined' to make an all-libraries file.
      Example:
        -L combined,SL_FOS to make a combined set and a fosmid ends set.
      Defaults to combined, plus all available libs for the given organism.

   -O common_name,common_name
      list of organisms to include in the sequence release.  Case
      insensitive.
      Default 'tomato'

EOU
}
sub HELP_MESSAGE {usage()}

sub cmd_tee(@) {
    print "running: ",@_,"\n";
    return @_;
}


our %opt;
getopts('vL:O:',\%opt) or usage();
$opt{O} ||= 'Tomato';
my $org_sql = join ',', map qq|'$_'|, map lc, split /,/,$opt{O}; #< sql-quote the org names

if($ARGV[0]) {
  chdir $ARGV[0]
    or die "could not chdir to '$ARGV[0]': $!\n";
}

my $dbh = CXGN::DB::Connection->new({ dbargs => {AutoCommit => 1}});
my $accession_ids = $dbh->selectcol_arrayref(<<EOQ);
select a.accession_id
from sgn.accession a
join sgn.organism o using(organism_id)
join sgn.common_name cn using(common_name_id)
where lower(cn.common_name) IN($org_sql)
EOQ

#warn "$opt{O}\n", Dumper $accession_ids;
#die;

my @shortnames = $opt{L}
    ? split /,/,$opt{L}
    : ('combined',map {$_->shortname} map {CXGN::Genomic::Library->search( accession_id => $_)} @$accession_ids);

#now validate all the shortnames
foreach my $sn (@shortnames) { 
  if( lc($sn) eq 'combined' ) {
    $sn = '';
  } else {
    my ($lib) = map CXGN::Genomic::Library->search( shortname => $sn, accession_id => $_ ), @$accession_ids;
    $lib or die "no library found for organisms $opt{O} with shortname '$sn'\n";
  }
}

$dbh->disconnect(42);

foreach my $libarg (@shortnames) {
  my $lib_fn = $libarg || 'combined';
  $libarg &&= "-L $libarg";

  system cmd_tee "genomic_get_seqs.pl -O $opt{O} -q $libarg -T -F  bacends_${lib_fn}_raw";
  system cmd_tee "genomic_get_seqs.pl -O $opt{O} -q $libarg bacends_${lib_fn}_screened_and_trimmed";
  make_blast_db_tar("bacends_${lib_fn}_screened_and_trimmed.seq");
  system cmd_tee "gzip -f bacends_${lib_fn}_raw.seq bacends_${lib_fn}_raw.qual";
  system cmd_tee "gzip -f bacends_${lib_fn}_screened_and_trimmed.seq bacends_${lib_fn}_screened_and_trimmed.qual";
}


############# SUBROUTINES ###########

sub make_blast_db_tar {
  my ($seqfile) = @_;

  my $bn = basename($seqfile,'.seq');

  mkdir $bn;
  system "formatdb -p F -l /dev/null -i $seqfile -n $bn/$bn";
  system "tar czf blastdb_$bn.tar.gz $bn";
  system "rm -rf $bn";
}

