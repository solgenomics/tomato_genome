#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use Data::Dumper;

use File::Copy;

use CXGN::TomatoGenome::BACSubmission;
use CXGN::TomatoGenome::BACPublish;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script -s sp_organization_id  [other options] <tarball>

  Retrofit a bunch of submission tarballs with sequencer information,
  based on the given sp_organization_id.

  Republishes the given file using bsub_publish_bac_seq.pl

  Options:

    -u upload_account
       set the organization upload account name

    -s organization_shortname
       set organization shortname of the organization to set as sequencer

    -d <dir>
       directory to publish to, passed through to bsub_publish_bac_seq.pl
       if present

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('s:u:',\%opt) or usage();

usage() unless @ARGV && ( $opt{s} || $opt{u} );

foreach my $subfile (@ARGV) {

  my $submission = CXGN::TomatoGenome::BACSubmission->open($subfile);

  #set the sequencer info and check that it took
  my %info = $submission->sequencing_info(
					  $opt{s} ? (org_shortname => $opt{s}) : (),
					  $opt{u} ? (org_upload_account_name => $opt{u}) : (),
					 );
  my @possible_orgs = $submission->sequencing_possible_orgs;
  unless( @possible_orgs ) {
    die "can't find any organizations matching ".Dumper(\%info);
  }

  #now republish the submission tarball
  eval {
    CXGN::Tools::Run->run( 'bsub_publish_bac_seq.pl',
			   '-f',
			   '-G',
			   $opt{d} ? ( -d => $opt{d} ) : (),
			   $submission->new_tarfile
			 );
  };
  if( $EVAL_ERROR ) {
    warn "$subfile failed:\n$EVAL_ERROR\n";
  }
}


