#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Temp qw/tempfile/;
use File::Spec;

use CXGN::DB::Connection;

use CXGN::Tools::Script qw/lock_script unlock_script dprint debugging/;
use CXGN::Tools::Run;

use CXGN::Publish qw/publish/;
use CXGN::TomatoGenome::BACPublish qw/aggregate_filename glob_pattern/;
use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident/;

use CXGN::TomatoGenome::Config;

#use Data::Dumper;

#########  DEFAULTS ##########

our $cfg = CXGN::TomatoGenome::Config->load_locked;
our $pubdir = File::Spec->catdir( @{$cfg}{'ftpsite_root','bac_publish_subdir'} );

######### /DEFAULTS ##########

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options]

  Script to create aggregate files of the sequences and annotation
  results in the BAC publishing repository.  Uses CXGN::Publish to
  publish the aggregate files.  This script is meant to be run by
  cron.

  Options:

  -g go. don't need this if you gave other options.

  -d <dir>

     Publishing directory location.  Defaults to
     $pubdir

  -x   dry run, just print what you _would_ do

  -v   be verbose about what you're doing

EOU
}

lock_script() or die "Don't run more than one instance of $FindBin::Script at once\n";

#parse and validate command-line options
our %opt;
getopts('gxvd:',\%opt) or usage();
%opt or usage;

$pubdir = $opt{d} if $opt{d};
-d $pubdir or die "Publishing directory '$pubdir' does not exist\n";
-w $pubdir or die "Publishing directory '$pubdir' is not writable\n";

$CXGN::Publish::dry_run   = 1 if $opt{x};
$CXGN::Publish::print_ops = 1 if $opt{v};

my @pub; #< array to hold all the publish operations we need to do

#make the aggregate files for all bacs, all finished bacs,
#and per-chromosome for the finished bacs on each chrom
foreach my $set ( 'all', 'finished', map{"chr${_}_finished"} (0..12)  )  {
  foreach my $tagtype (qw/seqs rm_seqs/) {
    push @pub, cat_aggregate($set.'_'.$tagtype);
  }
  #now aggregate the gff3 files
  push @pub, gff3_aggregate($set.'_gff3');

  my $agg_fn = aggregate_filename($set.'_accs',$pubdir);
  my $glob   = glob_pattern($set.'_accs',$pubdir);
  if($agg_fn and $glob and my $accsfile = extract_accessions_from_seqfiles($glob)) {
    push @pub,['cp',$accsfile,$agg_fn];
  }
}

#now do the publishing operations
publish(@pub);

#allow other copies of this script to be run now
unlock_script();

exit;


########### SUBROUTINES ############

# open each of the gff3 files matching the pattern given, aggregate
# them into a new gff3 file, and return a pub command to copy that
# file into the proper destination
sub gff3_aggregate {
  my $tag = shift;
  my $dest = aggregate_filename($tag,$pubdir);
  my $pat = glob_pattern($tag,$pubdir);

  my (undef,$tempfile) = tempfile(File::Spec->catfile(File::Spec->tmpdir,'bsub-aggregate-annotations-gff3-XXXXXXXX'), UNLINK => 1);

  my $files_temp = files_tempfile($pat);

  my $agger = CXGN::Tools::Run->run('gff3_reformat.pl',
				    '-U',
				    -e => '/[\._]CDS_|_match_part_|_exon_/i',
				    -f => $files_temp,
				    { out_file => $tempfile }
				   );
  return ['cp',$tempfile,$dest];
}

# cat the pattern given into a temp file, generate a pub command to
# copy that temp file into the given destination
sub cat_aggregate {
  my $tag = shift;
  my $dest = aggregate_filename($tag,$pubdir);
  my $pat = glob_pattern($tag,$pubdir);

  my (undef,$tempfile) = tempfile(File::Spec->catfile(File::Spec->tmpdir,'bsub-aggregate-annotations-cat-XXXXXXXX'), UNLINK => 1);

#  warn "catting $pat to $tempfile, copy to $dest\n";

  my $file_list = files_tempfile( $pat );
  my $cmd = "xargs -a $file_list cat > $tempfile";
  system $cmd;
  if($CHILD_ERROR) {
    warn "command failed ($CHILD_ERROR,$!): $cmd\n";
    return ();
  }
  return ['cp',$tempfile,$dest];
}

sub files_tempfile {
  my ($pat) = @_;

  my ($files_temp_fh,$files_temp) = tempfile(File::Spec->catfile(File::Spec->tmpdir,'bsub-aggregate-annotations-filelist-XXXXXXXX'), UNLINK => 1);

  my @files = glob($pat);

  $files_temp_fh->print("$_\n") foreach @files;

  close $files_temp_fh;

  return $files_temp;
}

# given a glob pattern specifying a bunch of sequence files, extract
# the SGN idents and genbank accessions from their deflines and return
# a tab-delimited tempfile listing them as SGN_ID\tACC SGN_ID\tACC ...
sub extract_accessions_from_seqfiles {
  my ($seqfiles_pat) = @_;

  my @seqfiles = glob $seqfiles_pat
    or return;

  my ($t_fh,$tempfile) = tempfile(UNLINK => 1);

  foreach my $seqfile ( @seqfiles ) {
    open my $deflines, "grep '>' $seqfile |"
      or die "could not run grep on file '$seqfile': $!";
    while(my $dl = <$deflines>) {
      my ($sgnid,$acc) = $dl =~ /^\s*>\s*(\S+)\s+(\S+)/;
      my $p = parse_clone_ident($sgnid,'versioned_bac_seq')
	or die "sequence name '$sgnid' is not a valid versioned bac seq name in file '$seqfile'";
      if($acc =~ /^[A-Z][A-Z_]{1,3}\d+\.\d+$/) {
	print $t_fh "$sgnid\t$acc\n";
      }
    }
  }
  return $tempfile;
}
