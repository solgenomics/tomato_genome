#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use File::Basename;

use CXGN::TomatoGenome::BACSubmission;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script  tar_file tar_file ...

  Go through all the sequence tarballs given on the command line,
  figure out their HTGS phase number (1, 2, or 3), and load it
  as a featureprop in our chado database.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();

@ARGV or usage;

-r or die "cannot open file $_\n" foreach @ARGV;


our $dbh = CXGN::DB::Connection->new({dbargs => {AutoCommit => 0}});

my $cvt_id = get_phase_cvterm_id()
  or die "got no cvt_id!";

my $q_insert_featureprop = $dbh->prepare(<<EOQ);
insert into featureprop (feature_id,type_id,value)
values ( (select feature_id from feature where uniquename=?),
         $cvt_id,
         ?
       )
EOQ

my $q_check_featureprop = $dbh->prepare(<<EOQ);
select fp.value
from featureprop fp
join feature f using(feature_id)
where f.name = ?
  and fp.type_id = $cvt_id
EOQ

foreach my $tarfile (@ARGV) {
  eval {
    my $sub = CXGN::TomatoGenome::BACSubmission->open($tarfile);
    my $bn = basename($tarfile);
    if ( my $e = $sub->validation_text ) {
      die $e;
    }

    my $seqname = $sub->sequence_identifier;

    #check if this sequence already has a phase on its feature
    $q_check_featureprop->execute($seqname);
    my ($db_phase) = $q_check_featureprop->fetchrow_array;

    my $sub_phase = $sub->htgs_phase;
    die "invalid sub phase $sub_phase" unless $sub_phase == 1 || $sub_phase == 2 || $sub_phase == 3;

    if ($db_phase && $db_phase == $sub_phase) {
      die "db already has phase $db_phase. skipping\n";
    }

    $db_phase && !($db_phase == $sub_phase)
      and die "db phase $db_phase != submission phase $sub_phase";

    warn "$seqname: setting htgs_phase $sub_phase\n";
    $q_insert_featureprop->execute($seqname,$sub_phase);
  };
  if($EVAL_ERROR) {
    warn "skipping $tarfile:\n$EVAL_ERROR";
    $dbh->rollback;
  } else {
    $dbh->commit;
  }
}

sub select_phase_cvterm_id {
  my ($htgs_phase_cvterm_id) = $dbh->selectrow_array(<<EOQ); #<guaranteed by unique constraint to return either 0 or 1 rows
select cvterm_id
from cvterm ct
join dbxref using(dbxref_id)
where ct.name = 'htgs_phase'
  and dbxref.accession='autocreated:htgs_phase'
  and ct.is_obsolete = 0
order by cvterm_id desc
limit 1
EOQ

  return $htgs_phase_cvterm_id;
}

sub get_phase_cvterm_id {

  my $cvtid = select_phase_cvterm_id();
  return $cvtid if $cvtid;

#  die "no htgs_phase cvterm id found";

  #otherwise, we have to create it....

  warn "creating new htgs_phase featureprop type\n";

  my $propname = 'htgs_phase';
  my $dbxname = "autocreated:$propname";

   $dbh->do(<<EOQ,undef,$dbxname);
 insert into dbxref (db_id, accession)
 values ( (select db_id from db where name='null' and description like '%fake %'),
          ?
        )
EOQ
  my ($dbx_id) = $dbh->selectrow_array("select dbxref_id from dbxref where accession = ?",undef,$dbxname);
  die "no dbxref_id found" unless $dbx_id;
  warn "got dbxref_id $dbx_id\n";

  $dbh->do(<<EOQ,undef,$propname,$dbx_id);
insert into cvterm (cv_id,name,dbxref_id)
values ( (select cv_id from cv where name = 'feature_property'),
          ?,
          ?
       )
EOQ

  return select_phase_cvterm_id() || die "still no cvterm id!";
}
