#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Spec;

use Data::Dumper;

use CXGN::Genomic::Clone;
use CXGN::Genomic::Chromat;

use CXGN::Tools::List qw/balanced_split collate min all/;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] -f fpc_file  bac_names_file

  Given a list of BAC names, filter them for only BACs with two good
  end sequences, not yet assigned for IL mapping, not yet sequenced,
  and with BAC ends that have not hit any sequenced BAC, then
  partition them equally among the tomato genome sequencing projects
  with IDs specified with -P.

  Produces files IL_bacs.<country>.txt in the output dir specified by
  -o.

  Options:

  -p <id>
      set the sp_person_id to use for loading the assignments,
      if -L was specified.  Default 290 (Rob Buels).

  -L  after generating the assignments, load them into the DB

  -s <num>
      set max number of BACs to assign per project.  Default 200.

  -o <dir>
      set the directory to put the IL_bacs.*.txt files.  Default '.'

  -P <num>,<num>,<num>
      list the project ids to divide the BACs among.
      Default 1,2,3,4,5,6,7,8,9,12

  -f <file>
      GFF file containing FPC map results to use.

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('Ls:o:p:P:f:',\%opt) or usage();
@ARGV or usage;
$opt{s} ||= 200;
$opt{p} ||= 290;
$opt{f} && -f $opt{f} or usage('must give a gff file of fpc data with -f');

our @projects = $opt{P} ? split /,/,$opt{P} : (1..9,12);
for (@projects) {
  die "invalid project $_\n" unless $_ >0 && $_ <=12;
}
$opt{o} ||= '.';
-d $opt{o} && -w $opt{o}
  or die "please use -o to set a writable directory\n";

our %stats; #< keep some global counters for statistics
our ($nr_bacs_file) = @ARGV;

#load indexed list of non-repetitive bacs
our %non_repetitive_bacs = do{
  open my $nrbacs,$nr_bacs_file
    or die "could not open non-repetitive bacs list $ARGV[0]: $!";
  map {
    chomp;
    $_ => 1
  } <$nrbacs>
};


our @picked_bacs; #< array of all BACs we'll be assigning
sub enough_bacs {
  @picked_bacs >= @projects*$opt{s}
}


my %fpc_contigs;
my @fpc_singletons;
# go through all the bigger-than-one FPC contigs in the Sanger set,
# and pick one BAC from each one that consists of all good BACs
open my $fpc_gff,$opt{f}
  or die "could not open $opt{f} for reading: $!";
my $s = 0;
while(my $line = <$fpc_gff>) {
  my @f = split /\t/,$line;
  my $ctg = $f[0];
  my ($bac) = split /;/,$f[8];
  next unless $bac =~ /^BAC "([^"]+)"/;
  $bac = $1;

  if($ctg eq 'ctg0') {
    push @fpc_singletons,$bac;
  } else {
    $fpc_contigs{$ctg} ||= [];
    push @{$fpc_contigs{$ctg}},$bac;
  }
}

#go through all the contigs now
while(my ($ctg,$bacs) = each %fpc_contigs) {
  my @clones = map {$non_repetitive_bacs{$_} && lookup_clone($_)} @$bacs;

  next unless all @clones;
  next unless all map {is_good_bac($_)} @clones;

  push @picked_bacs,$clones[0];
  print "picked ".$clones[0]->clone_name." from $ctg\n";
  last if enough_bacs;
}

#pick up the fpc singletons if we don't have enough
unless( enough_bacs ) {
  foreach my $bac (@fpc_singletons) {
    my $clone = lookup_clone($bac)
      or next;
    next unless is_good_bac($clone);
    push @picked_bacs,$clone;
    last if enough_bacs;
  }
}

#now get ecoRI bacs if we still don't have enough
unless( enough_bacs ) {
  open my $nrbacs,$nr_bacs_file
    or die "cannot open $nr_bacs_file: $!";
  while(my $bacname = <$nrbacs>) {
    chomp $bacname;
    next unless $bacname =~ /EcoRI/;
    my $clone = lookup_clone($bacname)
      or next;
    push @picked_bacs,$clone if is_good_bac($clone);
    last if enough_bacs;
  }
}

#now divide the BACs among the projects
my %assignments = collate \@projects,[balanced_split(scalar(@projects),@picked_bacs)];
#clamp the jobs to the -s option
foreach my $v (values %assignments) {
  $v = [@{$v}[0..min($opt{s}-1,@$v-1)]];
}

#count the job size we've produced for each of the projects
foreach (@projects) {
  $stats{"project_job_size_$_"} = $assignments{$_} ? scalar @{$assignments{$_}} : 0;
}

#print Dumper \%assignments;
print Dumper(\%stats);

#now load the assignments into the DB
if($opt{L}) {
  CXGN::Genomic::Clone->db_Main->dbh_param(AutoCommit => 1);
  while ( my ($project_id,$bacs_set) = each %assignments ) {
    foreach my $clone (@$bacs_set) {
      print "assigning ".$clone->clone_name." to project $project_id\n";
      $clone->il_mapping_project_id($project_id,$opt{p});
    }
  }
  CXGN::Genomic::Clone->db_Main->dbh_param(AutoCommit => 0);
}

#print out the assignments
while ( my ($project_id,$bacs_set) = each %assignments ) {
  my @projmap = qw/dummy USA Korea China UK India Netherlands France Japan Spain USA USA Italy/;
  my $filename = File::Spec->catfile($opt{o},"IL_bacs.$projmap[$project_id].txt");
  open my $projfile, ">$filename"
    or die "could not open $filename for writing: $!";

  foreach my $clone (@$bacs_set) {
    print $projfile $clone->clone_name, "\n";
  }
}


sub lookup_clone {
  my ($bacname) = @_;
  my $clone = CXGN::Genomic::Clone->retrieve_from_clone_name($bacname)
    or do{ print "no clone found for '$bacname', skipping\n";
	   $stats{'name not found'}++;
	   return 0;
	 };
  return $clone;
}

sub is_good_bac {
  my ($clone) = @_;
  my $bacname = $clone->clone_name;

  unless($non_repetitive_bacs{$bacname}) {
    print "$bacname is not in the non-repetitive BACs set\n";
    $stats{'not in non-repetitive'}++;
    return 0;
  }

  #skip this BAC unless it has at least two good reads, and at least
  #one from each end
  my @chromats = grep {
    my @g = grep {
      0 == keys %{$_->flags}
    } $_->gss_objects;
    @g >= 1 #< must have at least 1 non-flagged gss
  } $clone->chromat_objects;

  unless(@chromats >= 2) {
    print "only ".scalar(@chromats)." reads for $bacname, skipping\n";
    $stats{'not enough reads'}++;
    return 0
  }
  unless(grep {$_->primer eq 'SP6'} @chromats) {
    print "no good SP6 read for $bacname, skipping\n";
    $stats{'no sp6'}++;
    return 0
  }
  unless(grep {$_->primer eq 'T7'} @chromats) {
    print "no good T7 read for $bacname, skipping\n";
    $stats{'no t7'}++;
    return 0
  }

  #skip this BAC if it's already assigned to a project in the DB
  if(my $pid = $clone->il_mapping_project_id) {
    print "$bacname is already assigned to project $pid, skipping\n";
    $stats{'already assigned'}++;
    return 0;
  }

  #skip this BAC if it's already being sequenced
  if(my $chr = $clone->chromosome_num) {
    print "$bacname is being sequenced on chromosome $chr\n";
    $stats{"sequencing on chr $chr"}++;
    return 0
  }

  #skip this BAC if any of its good bac ends hit already-sequenced
  #BACs
  my @aligned_bac_ends = grep {
    my $chromat_name = $_->clone_read_external_identifier;

    #look up hits on this one in the chado DB
    my ($cnt) = $_->db_Main->selectrow_array(<<EOQ,undef,$chromat_name);
select count(*)
from public.feature_synonym fs
join public.synonym s using(synonym_id)
where s.name = ?
EOQ
    $cnt
  } @chromats;

  if(@aligned_bac_ends) {
    print "$bacname has one or more ends aligned to other BACs, skipping\n";
    $stats{'aligned bac ends'}++;
    return 0;
  }

  #check that it has no fish data
  my ($fish) = $clone->db_Main->selectrow_array('select count(*) from sgn.fish_result where clone_id=?',undef,$clone->clone_id);
  if($fish) {
    print "$bacname has FISH data, skipping.\n";
    $stats{'has fish'}++;
    return 0;
  }

  #check that it has no overgo data
  my ($overgo) = $clone->db_Main->selectrow_array('select count(*) from physical.bac_associations where bac_id=?',undef,$clone->clone_id);
  if($overgo) {
    print "$bacname has overgo data, skipping.\n";
    $stats{'has overgo'}++;
    return 0;
  }

  print "$bacname is good\n";

  return 1;
}
