#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;
use File::Basename;

#use Data::Dumper;

use CXGN::Tools::List qw/all/;
use CXGN::TomatoGenome::BACSubmission;
use CXGN::TomatoGenome::BACPublish qw/publishing_locations/;
use CXGN::Publish qw/move_or_print/;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  my $analyses_list_string = join ', ', map {$_->analysis_name} CXGN::TomatoGenome::BACSubmission->analyses;
  my $params_list_string = join "\n", map { my %params = $_->list_params;
					    map { "* $_\n  $params{$_}" } keys %params
					  } CXGN::TomatoGenome::BACSubmission->analyses;

  $params_list_string =~ s/\n/\n     /g; #indent by 5 spaces
  my $t = File::Spec->tmpdir;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] <submission tarball> <submission tarball>

  Script to analyze a BAC submission using a subset of the analyses in
  the BAC pipeline, placing the analysis result files in the directory
  you specify with the -d option.

  Options:

  -a <analysis_name,analysis_name,...>
     Comma-separated list of analysis names to run on these submissions.
     Right now, available analyses are:
     $analyses_list_string

  -d  directory in which to put the result files

  -o tag=value,tag=value,tag=value,...
     List of parameters for the analyses you want to run.
     Params are:
     $params_list_string

  -t <dir>
     set base dir for temp files, default $t
EOU
}

### parse and validate arguments
our %opt;
getopts('a:d:o:t:',\%opt) or usage();

if( $opt{t} ) {
    $opt{t} = File::Spec->rel2abs( $opt{t} );
    CXGN::Tools::Wget->temp_root_dir( $opt{t} );
    CXGN::TomatoGenome::BACSubmission->tempdir( $opt{t} );
}

my @submission_tarballs = @ARGV
  or usage('No BAC submission tarballs specified');
$opt{a} or usage('Must specify at least one analysis name');
my @analysis_names = split /,/,$opt{a};
$opt{d} or usage('No -d option specified');
-w $opt{d} or die "Directory '$opt{d}' not found or not writable\n";
my %analysis_opts = map { my ($t,$v) = split /=/,$_;
			  $t && $v or die "Invalid tag/value pair '$_'";
			  $t,$v
			} split /,/,($opt{o}||'');

#check that all analysis names they requested actually exist
foreach my $aname (@analysis_names) {
    CXGN::TomatoGenome::BACSubmission->get_analysis($aname)
          or die "Unknown analysis name '$aname'.  Available analyses are ".join(',',CXGN::TomatoGenome::BACSubmission->list_analyses)."\n";
}

#check that the destination directory exists
#check that all the submission tarballs are real files
-r or die "Cannot open file '$_' for reading\n" for @submission_tarballs;

foreach my $subfile (@submission_tarballs) {
  my $sub = CXGN::TomatoGenome::BACSubmission->open_stripped($subfile);
  foreach my $aname (@analysis_names) {
    my $pub_locations = publishing_locations($opt{d},$sub->bac_name,$sub->is_finished);
    my $analysis_locations = $pub_locations->{"annot_$aname"}
      or die "No publishing locations found for analysis $aname!  This analysis needs to be added to the CXGN::TomatoGenome::BACPublish::publishing_locations function.\n";
    my @resultfiles = $sub->analyze_with($aname,\%analysis_opts);
    @$analysis_locations == @resultfiles
      or die "List of publishing locations does not match list of result files!\n";

    #flatten the publishing locations to just be the filenames
    @$analysis_locations = map { my $bn = basename($_); File::Spec->catfile($opt{d},$bn) } @$analysis_locations;

    while(my $rfile = shift @resultfiles) {
      my $dest = shift @$analysis_locations;
      move_or_print($rfile,$dest);
    }
  }
}


