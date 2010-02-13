#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use File::Spec;

#use Data::Dumper;
use CXGN::Tools::Run;

our $ts_dir = '/data/prod/ftpsite/tomato_genome/bacs/training_set';

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options]

  Generate a training set from manually annotated gamexml files.

  Options:

  -g go.

  -d <dir>
    directory in which to find the training set directory hierarchy,
    default $ts_dir

  -r if set, remove redundant sequences in the gene model fasta files

  -S skip splitting the annotated BACs into gene models (using the cluster)

EOU
}

our %opt;
getopts('gd:rS',\%opt) or usage();
$opt{g} or usage();
$ts_dir = $opt{d} if $opt{d};
-d $ts_dir or die "$ts_dir is not a directory\n";
our $gene_models = "$ts_dir/gene_models";
our $annotated_bacs = "$ts_dir/annotated_BACs";

my @bac_xml_files = grep {-f} glob("$annotated_bacs/xml/*");

#extract the annotated bacs as gff 2 and 3 in the background
my @gff_jobs = map {
  CXGN::Tools::Run->run_async('gamexml_to_gff.pl',
				-v => $_,
				-o => "$annotated_bacs/gff/all.gff$_",
				'-C',
				@bac_xml_files,
			       );
} (2,3);

unless($opt{S}) {
#delete all the existing gene models
  CXGN::Tools::Run->run("rm -f $gene_models/xml/*");

  #check on our gff_jobs, in case they had some error
  $_->alive foreach @gff_jobs;

  #run gamexml_subregion on the cluster to extract the gene models
  my @split_jobs = map {
    CXGN::Tools::Run->run_cluster('gamexml_subregion.pl',
				  qw/-C -s/,
				  -o => "$gene_models/xml",
				  $_,
				  {
				   temp_base => '/data/shared/tmp', working_dir => '/data/shared/tmp' },
				 );
  } @bac_xml_files;

  #wait for the cluster jobs to finish
  sleep 1 while grep {$_->alive} @split_jobs;

  #clean up all the cluster jobs
  $_->cleanup foreach @split_jobs;
}

#check on our gff_jobs, in case they had some error
$_->alive foreach @gff_jobs;

my @fastadumps = map {
  my $spliced = $_;
  my $filename = $spliced ? 'spliced' : 'genomic';

  #convert the gene models to fasta
  CXGN::Tools::Run->run_async('gamexml_to_fasta.pl',
			      '-u',
			      $spliced ? '-s' : (),
			      -o => "$gene_models/fasta/$filename.seq",
			      grep {-f} glob("$gene_models/xml/*"),
			     );
} (0,1);

foreach (@fastadumps) {
  $_->wait;
  $_->cleanup;
}

my @nrs = map {
  my $spliced = $_;
  if($opt{r}) {
    my $filename = $spliced ? 'spliced' : 'genomic';
    CXGN::Tools::Run->run('fasta_nonredundant.pl',
			  -o => "$gene_models/fasta/${filename}_nonredundant.seq",
			  '-v',
			  "$gene_models/fasta/$filename.seq",
			 );
  } else {
    ()
  }
} (0,1);

foreach (@gff_jobs,@nrs) {
  $_->wait;
  $_->cleanup;
}
