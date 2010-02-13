#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Basename;
use File::Copy;

#use Data::Dumper;

use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident/;
use CXGN::Tools::Script qw/in_fh/;
use CXGN::TomatoGenome::BACSubmission;
use CXGN::TomatoGenome::BACPublish qw/publisher/;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] submission_file submission_file ...

  Set the GenBank accessions on a set of BAC submission tarballs.
  Alters the tarballs in-place (via temp files).  Takes the accessions
  to set either from a file specified with the -a option or from
  standard input.  By default, used CXGN::Publish when updating the
  file in-place, which saves the old version of the file and
  increments its file version number.

  Options:

  -a <filename>
     2-column tab-delimited file like:
     clone_name\taccession\n
     clone_name\taccession\n
     ...

  -f force overwriting existing GenBank accessions in the submissions

  -P do NOT use CXGN::Publish for updating the file, just modify the
     file in-place with no file versioning or backups

EOU
}

our %opt;
getopts('a:fP',\%opt) or usage();

@ARGV or usage;
-r or die "'$_' is not readable file\n" foreach @ARGV;

#parse and validate the genbank accessions
my %accessions;
my $acc_fh = in_fh($opt{a});
warn "reading accession assignments from STDIN.\n" if $acc_fh == \*STDIN;
my @acc_errors;
while(my $line = <$acc_fh>) {
  chomp $line;
  next unless $line =~ /\S/; #skip blank lines
  my ($clone_name,$accession) = split /\s+/,$line;
  parse_clone_ident($clone_name,'agi_bac_with_chrom')
    or push @acc_errors,"clone name '$clone_name' does not appear to be properly formatted";
  $accession !~ /\.\d+$/
    or push @acc_errors, "accession '$accession' is versioned.  Please use only unversioned identifiers.";
  $accession =~ /^[A-Z_]+\d+$/
    or push @acc_errors, "'$accession' doesn't look like a properly-formed genbank accession (e.g. AC171733)";
  $accessions{$clone_name} = $accession;
}
close $acc_fh;
#if there were errors in the accessions mapping, die and print them
#out
die map {"- $_\n"} @acc_errors if @acc_errors;

foreach my $bsub_file (@ARGV) {
  my $bn = basename($bsub_file,'.tar.gz');

  #skip this submission file unless it's one of the accessions we're
  #supposed to be setting
  unless(grep {index($bsub_file,$_) != -1} keys %accessions) {
    warn "skipping $bn, no accession given for it\n";
    next;
  }

  #open the submission, could take a long time to decompress it all
  my $submission = CXGN::TomatoGenome::BACSubmission->open($bsub_file)
    or die "could not open submission '$bsub_file'";

  #if it already has an accession, and we haven't been given the -f option,
  #then warn and skip it
  if(!$opt{f} and my $acc = $submission->genbank_accession) {
    warn "$bn already contains genbank accession '$acc', cowardly refusing to overwrite.  Use -f option to force.\n";
    next;
  }

  #now look up the accession we're supposed to be setting on it
  my $acc = $accessions{$submission->clone_object->clone_name_with_chromosome}
    or die "sanity check failed.  I thought I had a request to set an accession for the BAC in '$bn', but now I can't find it";

  #set the accession in the BAC submission, without the version number
  $submission->genbank_accession($acc,'do_not_append_version');

  #warn of any validation errors
  warn "WARNING, errors found after setting accession - ".$submission->validation_text if $submission->validation_errors;

  #make a new tarball and copy it over top of the original
  my $new_tar = $submission->new_tarfile;
  if($opt{P}) {
    move($new_tar,$bsub_file)
      or die "could not copy new tarball '$new_tar' => '$bsub_file': $!";
  } else {
    my $parsed = publisher->parse_versioned_filepath($bsub_file);
    my $pubdest = $parsed && $parsed->{full_unversioned} ? $parsed->{full_unversioned} : $bsub_file;
    publisher->publish([cp => $new_tar => $parsed->{full_unversioned}]);
#    warn join(',',cp => $new_tar => $parsed->{full_unversioned});
  }

  #clean up temp storage
  $submission->close;
}


