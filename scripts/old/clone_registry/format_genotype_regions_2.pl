#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use Data::Dumper;
use CXGN::DB::Connection;

use CXGN::Marker::Tools qw/clean_marker_name marker_name_to_ids/;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script file file

  One-off script for initial loading of phenome.genotype_region table
  with IL polymorphic fragments and mapping bins derived from them.

  Expects CSV files on stdin or as command-line arguments.  CSV files
  should be formatted as:

  "marker name","chromosome num","line def","line def",...

  where a line def is a string like "IL1-2 N", which would indicate
  that the given marker is the northernmost known marker on IL1-2.

  Prints SQL on standard output that inserts the necessary rows into
  phenome.genotype phenome.genotype_region using something like:

  Options:

  -p <num>
    sp_person_id to use, default 290, which is Rob Buels

  -m string
    short name sql LIKE pattern for the map to use when loading
    default '\%EXPEN%2000%'

  -e string
    short name sql LIKE pattern for the genotype_experiment to
    use when loading
    default '\%EXPEN%2000%'

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('p:m:e:',\%opt) or usage();
$opt{p} ||= 290;
$opt{m} ||= '%EXPEN%2000%';
$opt{e} ||= '%EXPEN%2000%';

my $dbh = CXGN::DB::Connection->new;

#look up the map_version_id we'll need
my $map_matches= $dbh->selectall_arrayref(<<EOS,undef,$opt{m});
select map_version_id
from sgn.map_version mv
join sgn.map m using(map_id)
where m.short_name like ?
  and mv.current_version = true
EOS

@$map_matches == 1
  or die scalar(@$map_matches)." have shortnames matching '%opt{m}'";

my $map_version_id = $map_matches->[0][0];

my @stored_north;
my @stored_south;


my ($experiment_id) = $dbh->selectrow_array(<<EOSQL,undef,$opt{e});
select genotype_experiment_id
from phenome.genotype_experiment
where experiment_name ilike ?
EOSQL
$experiment_id or die "could not find genotype_experiment matching -e '$opt{e}'\n";


my @rows;

-f $_ && -r $_ or die "'$_' not an openable file\n"foreach @ARGV;

#snarf the whole file into memory, looking up some stuff as we do so
my @errors;
while(<>) {
  chomp;
  next unless $_;

  my ($markername,$chr,@ildefs) = map {s/^"|"$//g; $_} split ',',$_;
  die "no chr column" unless $chr;

#  warn "got ".join(',',$chr,$ident,@mns)."\n";

  my @ids = marker_name_to_ids($dbh,$markername)
    or push @errors, "no marker '$markername' found";
  push @errors, "weird ids $markername -> ".join(',',@ids)
    if @ids > 1;


  #look up the lg_id for this chromosome on the given map
  my ($lg_id) = $dbh->selectrow_array(<<EOS,undef,$map_version_id,$chr);
select lg_id from sgn.linkage_group where map_version_id = ? and lg_name = ?
EOS
  push @errors, "cannot find lg_id for chr $chr on map_version_id $map_version_id\n"
    unless $lg_id;

#  push @rows,[$markername,$lg_id,@ildefs];
  push @rows,[$ids[0],$lg_id,@ildefs];
#  print "made row ".Dumper($rows[-1]);
}


die "errors in input file:\n",map "  - $_\n",@errors if @errors;

my %il_defs;
for(my $i=0;$i<@rows;$i++) {
  my ($mid,$lgid,@ildefs) = @{$rows[$i]};

  sub get_maybe(@) {
    my ($r,$ind) = @_;
    return '' unless $r;
    return $r->[$ind];
  }

  foreach my $ildef (grep $_,@ildefs) {
    if( my ($il,$ns) = $ildef =~ /^(\S+)\s+(NS|N|S)/ ) {

      if($ns eq 'N' || $ns eq 'NS') {
	if($il =~ /^IL/) {
	  my $insert_genotype = $dbh->prepare_cached(<<EOSQL);
insert into phenome.genotype
 ( genotype_experiment_id, individual_id )
values
 ( (select genotype_experiment_id from phenome.genotype_experiment
    where experiment_name ilike ? and experiment_name like 'IL mapping%'),
   (select individual_id from phenome.individual where name = ?)
 )
EOSQL
	  my $get_genotype = $dbh->prepare_cached(<<EOS);
select genotype_id
from phenome.genotype g
join phenome.individual i
  using(individual_id)
join phenome.genotype_experiment ge
  using(genotype_experiment_id)
where
      ge.experiment_name like 'IL mapping%'
  and ge.experiment_name ilike ?
  and i.name = ?
EOS
	  $insert_genotype->execute($opt{e},$il);
	  $get_genotype->execute($opt{e},$il);
	  my ($genotype_id_il) = $get_genotype->fetchrow_array;

	  $il_defs{$il} = ['inbred', 'b', $genotype_id_il, $rows[$i]->[1], '', $opt{p}, get_maybe($rows[$i-1],0),$rows[$i]->[0]];

	  my $ilh = $il;
	  $ilh =~ s/^IL/ILH/;
	  $insert_genotype->execute($opt{e},$ilh);
	  $get_genotype->execute($opt{e},$ilh);
	  my ($genotype_id_ilh) = $get_genotype->fetchrow_array;
	  $get_genotype->finish;
	  $genotype_id_ilh or die "cannot find genotype_id for ILH '$ilh'\n";
	  $il_defs{$ilh} = [ 'inbred', 'h', $genotype_id_ilh,  $rows[$i]->[1], '', $opt{p}, get_maybe($rows[$i-1],0),$rows[$i]->[0]];

	} else {
	  $il_defs{$il} = ['bin', 'b', '', $rows[$i]->[1], $il, $opt{p}, get_maybe($rows[$i-1],0),$rows[$i]->[0]];
	}
      }

      if ($ns eq 'S' || $ns eq 'NS') {
	#S end of a definition
	push @{$il_defs{$il}}, $rows[$i]->[0], get_maybe($rows[$i+1],0);
	my $ilh = $il;
	if( $ilh =~ s/^IL/ILH/ ) {
	  push @{$il_defs{$ilh}}, $rows[$i]->[0], get_maybe($rows[$i+1],0);
	}
      }
    }
  }
}


#insert the necessary genotype_region rows
$dbh->do(<<EOSQL);
COPY phenome.genotype_region (type,zygocity_code,genotype_id,lg_id,name,sp_person_id,marker_id_nn,marker_id_ns,marker_id_sn,marker_id_ss) from stdin
EOSQL
foreach my $ilname (sort keys %il_defs) {
  my $row = $il_defs{$ilname};
#  print join("\t",$ilname, @$row),"\n";
  $dbh->pg_putline( join("\t",map {length($_) ? $_ : '\N'} @$row)."\n" );
}
$dbh->pg_endcopy;


$dbh->commit;
