#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Find;
use File::Basename;

use CXGN::DB::Connection;
#CXGN::DB::Connection->verbose(0); #STFU

use CXGN::Metadata;
use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident/;
use CXGN::Genomic::Clone;
use CXGN::Genomic::Search::Clone;

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message .= "\n" if $message;
  die $message,<<EOU;
Usage:
  $FindBin::Script -d <publishing directory>

  Check that all sequenced BACs are actually marked as such in the
  SGN database.

  Prints discrepancies on stdout in CSV format.

  Options:

  -d    base directory where sequenced bacs are published

EOU
}

our %opt;
getopts('d:',\%opt) or usage();
$opt{d} && -d $opt{d} or usage('Invalid -d option.');

my %sequenced_bacs = map {
  my ($bacname,$dirname) = fileparse($_,qr|\..+$|)
    or return;
  my $parsed = parse_clone_ident($bacname,'agi_bac_with_chrom')
    or return;
  #        $parsed->{match} eq $bacname
# 	 or return;
  $parsed->{dirname} = $dirname;

  #return
  ($bacname => $parsed)
} glob($opt{d}.'/chr*/*finished/*.tar.gz');

my @discrepancies;

my @sequenced_bac_names = sort keys %sequenced_bacs;

foreach my $bacname (@sequenced_bac_names) {
  my $clone = CXGN::Genomic::Clone->retrieve_from_parsed_name($sequenced_bacs{$bacname})
    or die "Could not retrieve clone $bacname\n";
  $sequenced_bacs{$bacname}->{clone_object} = $clone;
  my $dirname = $sequenced_bacs{$bacname}->{dirname};
  my $truestatus = $dirname =~ m|unfinished| ? 'in_progress' : 'complete';
  my ($chr_from_dir) = $dirname =~ /chr(\d{2})/;
  $chr_from_dir += 0; #make sure it's numeric
  $sequenced_bacs{$bacname}->{true_status} = $truestatus;
   unless( $clone->sequencing_status eq $truestatus ) {
    push @discrepancies,
      [ $clone->clone_id,
	$clone->chromosome_num || $chr_from_dir,
	$clone->clone_name_with_chromosome || $clone->clone_name,
	$truestatus,
	$clone->sequencing_status,
	$dirname,
      ];
  }
#   else {
#     print "actually right!\n";
#   }
}

#now find all the BACs that are marked complete in the database and aren't
#found on the filesystem
my $search = CXGN::Genomic::Search::Clone->new;
my $query = $search->new_query;
$query->sequencing_status('=?','complete');
my $results = $search->do_search($query);
$results->autopage($query,$search);
while(my $clone = $results->next_result) {
  my $chrname = $clone->clone_name_with_chromosome;
#  print $clone->clone_name." = $chrname\n";
  unless( my $record = $sequenced_bacs{$chrname} ) {
    push @discrepancies,
      [
       $clone->clone_id,
       $clone->chromosome_num || '',
       $chrname,
       'not present',
       $clone->sequencing_status,
       '',
      ];
  }
#   else {
#     print "actually right!\n"
#   }
}

print join("\t",'ID','Chromosome','Name','On FTP','In Database','Dir'),"\n";
foreach my $row (sort {$a->[1] <=> $b->[1]} @discrepancies) {
  print join("\t",@$row),"\n";
}



