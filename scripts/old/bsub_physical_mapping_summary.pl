#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use CXGN::DB::Connection;
use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident assemble_clone_ident/;
use CXGN::TomatoGenome::BACPublish qw/ aggregate_filename publisher /;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script

  Dump a tab-delimited file to STDOUT summarizing physical mapping
  info for all tomato BACs.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();

my $dbh = CXGN::DB::Connection->new();

my $orgname = 'Tomato';
warn "finding '$orgname' BAC libraries...\n";
my $library_ids = $dbh->selectcol_arrayref(<<EOS,undef,$orgname);
select
   l.library_id
from genomic.library l
     join public.organism using(organism_id)
where common_name  = ?
EOS

$library_ids = join ',',@$library_ids;

my $main_query = <<EOSQL;
-- Overgo data
select
	c.clone_id,
	(l.shortname || lpad(cast(c.platenum as text),4,'0') || c.wellrow || lpad(cast(c.wellcol as text),2,'0')) as clone_name,
        'Overgo' as method,
        ov.alias as marker,
	ov.mapname as map,
	ov.chromosome as chromosome,
	cast(ov.position as text) as position
from genomic.clone c
join genomic.library l using(library_id)
join ( 	select oa.bac_id as clone_id, ma.alias, pm.marker_id, mm.lg_name as chromosome, map.short_name as mapname, mm.position as position
		from physical.overgo_associations as oa
		join physical.oa_plausibility oap ON (oa.overgo_assoc_id = oap.overgo_assoc_id)
		join physical.probe_markers pm using(overgo_probe_id)
		join sgn.marker_to_map mm on (mm.marker_id = pm.marker_id)
                join sgn.marker_alias ma on (ma.marker_id = mm.marker_id and ma.preferred = true)
		join sgn.map map on (map.map_id = mm.map_id)
		where mm.current_version = true
	          and oap.plausible > 0
          ) as ov
      on (ov.clone_id = c.clone_id)
where l.library_id IN($library_ids)
UNION
-- IL mapping data
select
	c.clone_id,
	(l.shortname || lpad(cast(c.platenum as text),4,'0') || c.wellrow || lpad(cast(c.wellcol as text),2,'0')),
        'IL mapping',
        'IL_mapping',
        'Zamir IL Bins',
	il.chromosome || '',
	gr.name
from genomic.clone c
join genomic.library l using(library_id)
left join sgn_people.clone_il_mapping_bin_log il on (il.is_current = true and c.clone_id = il.clone_id)
left join phenome.genotype_region gr using(genotype_region_id)
where l.library_id IN($library_ids)
UNION
-- FISH data
select
	c.clone_id,
	(l.shortname || lpad(cast(c.platenum as text),4,'0') || c.wellrow || lpad(cast(c.wellcol as text),2,'0')),
        'FISH' as method,
        'FISH' as marker,
        'FISH' as map,
	fish.chromo_num || '' as chromosome,
	fish.chromo_arm || ' ' || fish.dist_avg || ' stdev ' || fish.dist_stddev as position
from genomic.clone c
join genomic.library l using(library_id)
join ( SELECT clone_id,
		   avg (percent_from_centromere) as pct,
             	   stddev (percent_from_centromere) as pct_stddev,
             	   avg (percent_from_centromere)*arm_length as dist_avg,
             	   stddev (percent_from_centromere)*arm_length as dist_stddev,
             	   arm_length,
             	   chromo_num,
           	   chromo_arm
	     FROM sgn.fish_result
	     NATURAL JOIN sgn.fish_karyotype_constants
	     GROUP BY clone_id,arm_length, chromo_num, chromo_arm
          ) as fish
	  on ( fish.clone_id = c.clone_id )
where l.library_id IN($library_ids)
ORDER BY clone_id
EOSQL

warn "running main query...\n";
my $sth = $dbh->prepare($main_query);
$sth->execute();
warn "done. gathering additional data and dumping...\n";

print join "\t", map qq|"$_"|,
  ( 'Clone ID',
    'Clone Name',
    'Method',
    'Marker',
    'Map',
    'Chromosome',
    'Position',
  );
print "\n";

my $prev_clone_id = 0;
while(my $row = $sth->fetchrow_arrayref() ) {
  my ($clone_id,$clone_name,@data) = @$row;

  my @rows_to_print;
  unless( $clone_id == $prev_clone_id ) {
    push @rows_to_print,
      get_marker_blast( $clone_id, $clone_name),
      get_sanger_fpc( $clone_id, $clone_name );
  }
  $prev_clone_id = $clone_id;
  push @rows_to_print, $row;

  my $prev_line = '';
  foreach (@rows_to_print) {
    my (undef,undef,@data) = @$_;
    next unless defined $data[-1];
    my $line = join("\t", map $_||'',@$_)."\n";
    print $line unless $prev_line eq $line;
    $prev_line = $line;
  }
}

warn "done.\n";

exit 0;

sub get_sanger_fpc {
  my ($clone_id, $clone_name) = @_;

  our $sanger_fpc_idx ||= do {
    print STDERR "indexing Sanger Tomato FPC results gff...";
    my %idx;
    my $fpc_gff = '/data/local/cxgn/core/sgn/documents/gbrowse/databases/fpc/tomato_R12_dQ.gff';
    open my $f, $fpc_gff or die "$! opening '$fpc_gff'.  Did you move it?";
    while(my $l = <$f>) {
      next if $l =~ /^\s*#/;
      my @f = split /\t/,$l;
      my $ctg = $f[0];
      next if $ctg eq 'ctg0'; #< ignore this
      my ($bac) = $f[8] =~ /BAC "(\w+)"/
	or next;
      $bac = _normalize_bac_name($bac)
        or next;

      $idx{$bac} = ['Sanger FPC','FPC','Sanger_FPC',undef,$ctg];
    }
    close $f;
    print STDERR "done.\n";
    # now return the finished index
    \%idx
  };

  return [$clone_id, $clone_name, @{$sanger_fpc_idx->{$clone_name}}] if $sanger_fpc_idx->{$clone_name};
  return ();
}


sub get_marker_blast {
  my ($clone_id, $clone_name) = @_;

  our $marker_blast_idx ||= do {
    print STDERR "indexing marker BLAST matches...";
    my %idx;
    my $all_bacs_gff = publisher()->published_as( aggregate_filename('all_gff3') )
        or die "could not find BACs combined gff3 file.  be sure you're running on a machine with access to the tomato genome BACs repository\n";
    -r $all_bacs_gff->{fullpath} or die "could not read BACs combined gff3 file '$all_bacs_gff->{fullpath}'\n";
    open my $marker_blast, "grep GenomeThreader_SGN_markers $all_bacs_gff->{fullpath} |"
      or die "$! grepping '$all_bacs_gff'";
    my %seen_match;
    while( my $l = <$marker_blast> ) {
      next if $l =~ /^\s*#/;
      my ($bac,undef,$type,@other) = split /\t/,$l;
      my $attr = pop @other;
      my ($marker_seq) = $attr =~ /Target=(\S+)/
	or next;
      $bac = _normalize_bac_name($bac)
	or next;

      my ($mid) = $marker_seq =~ /sgn-m(\d+)/i;
      my $seen_key = $mid || $marker_seq;
      next if $seen_match{$seen_key};
      $seen_match{$seen_key} = 1;
      
      #now look up the chromosome and position of this marker on the expen-2000
      my @results = do {
	if( $mid ) {
	  # look up the info on this marker
	  my $s = $dbh->prepare_cached(<<EOSQL);
	    select 'BLAST_vs_marker_seqs', ma.alias, map.short_name as mapname, mm.lg_name as chromosome, mm.position as position
	    from sgn.marker_to_map mm
            join sgn.marker_alias ma on (ma.marker_id = mm.marker_id and ma.preferred = true)
	    join sgn.map map on (map.map_id = mm.map_id)
	    where mm.current_version = true and mm.marker_id = ?
EOSQL
          $s->execute( $mid );
          my $r = $s->fetchall_arrayref;
          $s->finish;
          if( $r ) {
            #warn "got marker '$r->[0]' for mid '$mid'";
            ( @$r )
          } else {
	    ()
          }
	} else {
	  ()
	}
      };

      if( $idx{$bac} ) {
	push @{$idx{$bac}}, @results;
      } else { 
	$idx{$bac} = \@results;
      }
    }
    print STDERR "done.\n";
    \%idx;
  };

  return $marker_blast_idx->{$clone_name}
      # add the clone_id and clone_name to the rows before returning them, if found
    ? (map {[$clone_id,$clone_name,@$_]} @{$marker_blast_idx->{$clone_name}})
    : ();
}


sub _normalize_bac_name {
  my ($bac) = @_;
  my $p = parse_clone_ident($bac)
    or return;
  $bac = eval { assemble_clone_ident( agi_bac => $p ) }
    or return;
  return $bac;
}
