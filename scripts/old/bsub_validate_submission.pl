#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Spec;
use File::Basename;

#use Data::Dumper;

use CXGN::TomatoGenome::BACSubmission;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  my $tempdir = File::Spec->tmpdir;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] bac_tarball bac_tarball ...

  Script to check the validity of a Tomato BAC submission tarball from a sequencing center.

  Options:

  -r   attempt an in-place repair of the tarball, saving the original
       with a .bak extension.

  -R   same as -r, but do not save a backup file

  -t <dir>
      Set the directory this script uses for temporarily decompressing
      submission files. Defaults to $tempdir.
EOU
}

our %opt;
getopts('rRt:',\%opt) or usage();

#get any temporary directory the user specified
!$opt{t} || -d $opt{t} && -w $opt{t}
  or die "$opt{t} is not a writable directory.\n";
CXGN::TomatoGenome::BACSubmission->tempdir($opt{t}) if $opt{t};

my @submission_files = @ARGV or usage;
-r or die "File '$_' not found.\n" foreach @submission_files;

foreach my $subfile (@submission_files) {
  my ($basename) = fileparse($subfile);
  my $warnings;
  my $sub;
  eval {
    $sub = CXGN::TomatoGenome::BACSubmission->open($subfile);
    if(my @errors = $sub->validation_errors ) {
      if($opt{r} || $opt{R}) {
	warn "attempting repair of $subfile\n";
	$sub->repair;
	@errors = $sub->validation_errors;
	unless(@errors) {
	  warn "success!\n";
	  unless($opt{R}) {
	    system('cp',$subfile,"$subfile.bak");
	    die 'error backing up current tarfile during repair' if $CHILD_ERROR;
	  }
	  system('cp',$sub->new_tarfile,$subfile);
	  die 'error copying repaired tarfile into position' if $CHILD_ERROR;
	} else {
	  warn "repair failed.\n";
	}
      }
      die map {"  - ERROR: ".$sub->error_string($_)."\n"} @errors if @errors;
    }
    warn "WARNING: $basename appears unfinished (not HTGS phase 3).\n" unless $sub->is_finished;
  }; if( $EVAL_ERROR ) {
    print "$basename failed:\n$EVAL_ERROR\n";
  } else {
    print "$basename passed.\n";
    print $warnings if $warnings;
    print "\n";
  }
  $sub->close if $sub;
}


