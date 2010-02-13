#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use Bio::DB::GFF;

use CXGN::DB::Connection;

use CXGN::DB::GFF::Versioned;

use CXGN::TomatoGenome::BACPublish qw/ aggregate_filename publisher /;
use CXGN::TomatoGenome::Config;
use CXGN::Tools::Script qw/lock_script unlock_script/;
use CXGN::IndexedLog;

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script

  This script makes a database named for this date, loads it with the
  latest gbrowse annotations, and sets its permissions.  Also deletes
  older gbrowse annotation snapshots.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();

lock_script() or die "Don't run more than one $FindBin::Script at once!\n";

my $dbh = CXGN::DB::Connection->new;

my ($seqs,$gff3) =
  map { my $p = publisher()->published_as($_) or die "no version of $_ seems to be published!"; $p->{fullpath} }
  map aggregate_filename($_),
  'all_seqs', 'all_gff3';

#have we loaded these already?
my $log = CXGN::IndexedLog->open( DB => $dbh, CXGN::TomatoGenome::Config->load_locked->{'bac_processing_log'} );
my $log_string = "RELOADED_GBROWSE_WITH $gff3";

unless( $log->lookup(content => $log_string) ) {

  my $bdb = CXGN::DB::GFF::Versioned->new( -db => 'cxgn_conf:bacs_bio_db_gff_dbname', -user => '' );

  $bdb->load_new($seqs,$gff3);

  $log->append($log_string);
}

unlock_script();
exit;

