#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use DBI qw/:sql_types/;
use CXGN::DB::Connection;
use CXGN::Tools::Text qw/list_to_string/;
#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script  library_shortname  plate_size  begin_plate-end_plate

  Example:  $FindBin::Script SL_FOS 384 1-144

  This script inserts records into the genomic.clone table as necessary to ensure that
  rows are present for the given library

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();

@ARGV == 3 or usage;

my ($lib_shortname,$plate_size,$plate_range) = @ARGV;

#check that we know about this plate size
my %plate_schemes = (  384 => { rows    => ['A'..'P'],
				columns => [1..24],
			      },
		       96  => { rows    => ['A'..'H'],
				columns => [1..12],
			      },
		    );
$plate_schemes{$plate_size}
  or die "unknown plate size $plate_size, I only know about plate sizes ".list_to_string(sort keys %plate_schemes)."\n";

#check that the plate range is valid
my ($plate_s,$plate_e) = $plate_range =~ /^(\d+)-(\d+)$/;
$plate_s && $plate_e && $plate_s < $plate_e
  or die "invalid plate range '$plate_range', must be two integers in ascending order, like '144-168'\n";


our $dbh = CXGN::DB::Connection->new({dbargs => {AutoCommit => 0}});

#look up the library ID for that shortname
my ($library_id,$clone_type_id) = do {
  my $libs = $dbh->selectall_arrayref('select library_id,clone_type_id,shortname from genomic.library where lower(shortname)=lower(?)',undef,$lib_shortname);

  @$libs
    or die "no library found matching shortname '$lib_shortname'";

  @$libs > 1
    and die "multiple libraries found with shortname '$lib_shortname':\n"
      .( join '', map "  id: $_->[0]  shortname: $_->[2]\n", @$libs );

  $libs->[0][0], $libs->[0][1]
};


my $insert_clone = $dbh->prepare(<<EOQ);
insert into genomic.clone
         (  library_id,  clone_type_id, platenum, wellrow, wellcol )
  values ( $library_id, $clone_type_id,        ?,       ?,       ? )
EOQ
#index current lib clones by plate,row,col
my %curr_clones;
{
  my $c = $dbh->selectall_arrayref(<<EOQ);
select platenum || wellrow || wellcol
from genomic.clone
where library_id = $library_id
EOQ
  $curr_clones{$_->[0]} = 1 foreach @$c;
}

#now insert all the clones that need to be inserted
$dbh->do(<<EOQ);
COPY genomic.clone (  library_id,  clone_type_id, platenum, wellrow, wellcol )
FROM STDIN
EOQ

my $rows = $plate_schemes{$plate_size}{rows};
my $cols = $plate_schemes{$plate_size}{columns};
foreach my $plate ($plate_s..$plate_e) {
  foreach my $row (@$rows) {
    foreach my $col (@$cols) {
      unless($curr_clones{"$plate$row$col"}) {
	#warn "inserting $plate$row$col\n";
	$dbh->pg_putline("$library_id\t$clone_type_id\t$plate\t$row\t$col\n");
      }
#       else {
# 	warn "skipping $plate$row$col, already present\n";
#       }
    }
  }
}
$dbh->pg_endcopy;

$dbh->commit;
