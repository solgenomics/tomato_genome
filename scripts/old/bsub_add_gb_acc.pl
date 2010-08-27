#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use IO::Handle;
use File::Basename;

use CXGN::TomatoGenome::BACSubmission;
use CXGN::TomatoGenome::BACPublish qw/parse_filename bac_publish/;
use CXGN::Tools::List qw/all any/;
use CXGN::Tools::Run;
use CXGN::Genomic::CloneIdentifiers;
use CXGN::Tools::Identifiers qw/identifier_namespace clean_identifier/;

############### CONFIGURATION VARIABLES ##################

our $vhost = CXGN::TomatoGenome::Config->new;
our $country_uploads_path = $vhost->get_conf('country_uploads_path');
our $publish_path = File::Spec->catdir($vhost->get_conf('ftpsite_root'),$vhost->get_conf('bac_publish_subdir'));
our $logfile = File::Spec->catdir('/data/shared/tomato_genome/bacpipeline/processing_log');

##########################################################


#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script  -a list_file submission_file submission_file ...

  Given a bunch of tarballs on the command line and a list of GenBank
  accessions with -a, modify each submission in-place, adding the
  appropriate GenBank accession to it.

  Do not run this directly on published BAC submissions.  Copy them
  somewhere else, add the GenBank accessions, then republish them
  using bsub_publish_bac_submission.pl .

  Options:

    -a <file>
      accession list filename, containing a two-column list of BAC
      names and genbank accessions.  These must not contain version
      numbers.

    -B
      do not make .bak backup files of each submission that is
      modified

    -f
      force replacing accessions inside submissions that already have
      them

    -p
      instead of modifying it in-place, republish the modified files
      with bsub_publish_bac_submission.pl

    -d <dir>
       if republishing with -p, publish to this directory
       Default: $publish_path

    -l <file>
       if republishing with -p, log publishes to the given file
       Default: $logfile

    -x dry run, don't modify any files

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('a:Bfxp',\%opt) or usage();
$logfile = $opt{l} if $opt{l};
$publish_path = $opt{d} if $opt{d};
-d $publish_path or die "can't find publish dir '$publish_path'\n";

@ARGV && all map -f,@ARGV
  or usage('must provide some submission files');


### READ IN AND INDEX ALL THE ACCESSIONS ###
my %accessions;
open my $accs,$opt{a}
  or die "$! opening accessions from $opt{a}\n";
sub acc_err (@) {
  warn @_;
  our $accs_error = 1;
}
while(<$accs>) {
  my ($bacname,$accession) = split;
  my $line = $accs->input_line_number;
  identifier_namespace($bacname) eq 'bac'
    or acc_err "error at $opt{a} line $line - Invalid bac identifier '$bacname'\n";
  $bacname = clean_identifier($bacname,'bac');
  identifier_namespace($accession) eq 'genbank_accession' && $accession =~ /^[A-Z]{2}_?\d+$/
    or acc_err "error at $opt{a} line $line - Invalid genbank accession '$accession'\n";

  $accessions{$bacname} = { acc => $accession, used => 0, bac => $bacname };
}
our $accs_error and die "Errors in accessions file, aborting\n";

### NOW PROCESS EACH SUBMISSION
foreach my $subfile (@ARGV) {
  my $pfn = parse_filename($subfile); #< parsed filename
  our $bn = $pfn->{basename};

  sub subwarn(@) {
    my @p = @_;
    foreach(@p) {
      chomp;
      s/\n/\n$bn: /g;
    }
    warn map "$bn: $_\n",@p
  };

  unless( any map index($bn,$_) != -1, keys %accessions ) {
    subwarn "no accession given to add, skipping.\n";
    next;
  }

  my $sub = CXGN::TomatoGenome::BACSubmission->open($subfile);
  my $old_acc = $sub->genbank_accession;
  my $old_err_cnt = $sub->validation_errors;
  my $old_err_text = $sub->validation_text;
  my $acc_rec = $accessions{$sub->bac_name}
    or do { subwarn "no accession given to add to bac ".$sub->bac_name.", skipping";
	    next;
	  };
  $sub->genbank_accession($acc_rec->{acc});
  my $new_err_cnt = $sub->validation_errors;
  my $new_err_txt = $sub->validation_text;
  $acc_rec->{used} = 1;

  if($new_err_txt) {
    subwarn "still had errors after setting accession $acc_rec->{acc}";
    subwarn $new_err_txt;
  }

  if(!$opt{f} && $old_err_cnt < $new_err_cnt) {
    $old_err_text ||= "none\n";
    subwarn "old errors:\n$old_err_text";
    subwarn "new errors are more numerous than old errors, and -f not specified.  skipping.";
    subwarn "are you sure $acc_rec->{acc} is the right accession for $acc_rec->{bac}?";
    next;
  }

  unless($opt{x}) {
    #don't modify published submissions
    if($pfn->{file_version}) {
      if($opt{p}) {
	eval {
	  CXGN::Tools::Run->run('bsub_publish_bac_seq.pl',
				-d => $publish_path,
				'-f',
				-l => $logfile,
				$sub->new_tarfile,
			       );
	}; if($EVAL_ERROR) {
	  subwarn "publish failed";
	  subwarn $EVAL_ERROR;
	  next;
	}
      } else {
	subwarn "file appears to be versioned with CXGN::TomatoGenome::BACPublish, cowardly refusing to modify it in place.  maybe you want to use -p\n";
	next;
      }
    } else {
      $opt{B} or CXGN::Tools::Run->run( 'cp', $subfile, "$subfile.bak" );
      CXGN::Tools::Run->run( 'mv', $sub->new_tarfile, $subfile );
    }
  } else {
    subwarn "dry run, not writing change\n";
  }

  subwarn "finished setting $acc_rec->{acc}\n";
}
unless( all map $_->{used}, values %accessions ) {
  warn "WARNING: not all accessions in -a list were used:\n";
  warn map "  - $_->{bac} $_->{acc} not used\n", grep !$_->{used}, values %accessions;
}
