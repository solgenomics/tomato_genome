#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;
use CXGN::DB::Connection;

use CXGN::Marker::Tools qw/clean_marker_name marker_name_to_ids/;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script

  One-off script for initial loading of phenome.genotype_region table
  with IL polymorphic fragments and mapping bins derived from them.

  Options:

  -p <num>
    sp_person_id to use, default 290, which is Rob Buels

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('p:',\%opt) or usage();
$opt{p} ||= 290;

my $dbh = CXGN::DB::Connection->new;

#look up the map_version_id we'll need
my ($map_version_id) = $dbh->selectrow_array(<<EOS,undef,'%EXPEN%','%1992%');
select map_version_id
from sgn_bt.map_version mv
join sgn_bt.map m using(map_id)
where m.short_name like ? and m.short_name like ?
EOS

die 'could not find map_version_id for a map named like %EXPEN% %1992%'
  unless $map_version_id;

while(<>) {
  chomp;
  next unless $_;

  my ($chr,$ident,@mns) = map {s/^"|"$//g; $_} split "\t",$_;
  @mns = @mns[0..3];

#  warn "got ".join(',',$chr,$ident,@mns)."\n";

  die 'not enough marker names!' if @mns < 3;

  my @mids;
  foreach my $mn (@mns) {
    if($mn) {
      my @ids = marker_name_to_ids($dbh,$mn);
      die "weird ids $mn -> ".join(',',@ids)
	unless @ids == 1;
      push @mids,@ids;
    } else {
      push @mids,'';
    }
  }

  #look up the lg_id for this chromosome on the EXPEN 1992 map
  my ($lg_id) = $dbh->selectrow_array(<<EOS,undef,$map_version_id,$chr);
select lg_id from sgn.linkage_group where map_version_id = ? and lg_name = ?
EOS
  die "cannot find lg_id for chr $chr on map_version_id $map_version_id\n"
    unless $lg_id;

  sub row(@) {
    my @fields = @_;
    return (join "\t", map {qq|$_|} @fields)."\n";
  }

  #look up the genotype_id
  my $genotype_id_il;
  my $genotype_id_ilh;
  if( $ident  =~ /^IL/ ) {
    my $get_genotype = $dbh->prepare_cached(<<EOS);
select genotype_id
from phenome.genotype g
join phenome.individual i
  using(individual_id)
where experiment_name like 'IL-lines%'
  and name = ?
EOS
    $get_genotype->execute($ident);
    ($genotype_id_il) = $get_genotype->fetchrow_array;
    $get_genotype->finish;
    $genotype_id_il or die "cannot find genotype_id for IL '$ident'\n";
    print row 'inbred', 'b', $genotype_id_il,    $lg_id, '',    $opt{p}, @mids;

    my $ilh_ident = $ident;
    $ilh_ident =~ s/^IL/ILH/;
    $get_genotype->execute($ilh_ident);
    ($genotype_id_ilh) = $get_genotype->fetchrow_array;
    $get_genotype->finish;
    $genotype_id_ilh or die "cannot find genotype_id for ILH '$ilh_ident'\n";
    print row 'inbred', 'h', $genotype_id_ilh,  $lg_id, '',     $opt{p}, @mids;
  } else {
    print row 'bin',    'b', '',                $lg_id, $ident, $opt{p}, @mids;
  }
}
