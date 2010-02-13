#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use List::Util qw/max/;

use File::Basename;
use File::Spec;
use File::Temp qw/tempdir/;

use Bio::SeqIO;
use Bio::FeatureIO;

use CXGN::DB::Connection;

use CXGN::TomatoGenome::BACPublish qw/ agp_file contig_file /;
use CXGN::BioTools::AGP qw/ agp_to_seq agp_contig_seq /;
use CXGN::Tools::Script qw/lock_script unlock_script/;
use CXGN::Publish qw/ published_as publish /;
use CXGN::TomatoGenome::ChromosomeAssemblies qw/named_contigs contig_features/;
use CXGN::TomatoGenome::Config;

########### CONFIGURATION/DEFAULTS ################

my $cfg = CXGN::TomatoGenome::Config->load_locked;
our $agp_path = File::Spec->catdir( @{$cfg}{'ftpsite_root','agp_publish_subdir'} );
our $contigs_path = File::Spec->catdir(@{$cfg}{'ftpsite_root','contigs_publish_subdir'});

our @chromosome_nums = (0..12);

###################################################


sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options]

  Generates sequence sets according to the most recent published AGP
  files.

  Options:

    -a <dir>
      Directory to look in for AGP files.
      Default: $agp_path

    -d <dir>
      Directory to put contig files.
      Default: $contigs_path

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('a:d:',\%opt) or usage();

$contigs_path = $opt{d} if defined $opt{d};
$agp_path = $opt{a} if defined $opt{a};

#we will make our contig dir if not there
$CXGN::Publish::make_dirs = 1;
-d $contigs_path or warn "WARNING: contigs dir '$contigs_path' does not exist, will attempt to create it\n";
-d $agp_path or die "agp dir '$agp_path' does not exist\n";


lock_script() or die "please don't run more than one contig generation script at once.\n";

my $tempdir = tempdir(File::Spec->catdir(File::Spec->tmpdir, "$FindBin::Script-XXXXXX"), CLEANUP => 1);

my $all_ctg_tmp = File::Spec->catfile($tempdir,'all.contigs.fasta');
my $all_pm_tmp  = File::Spec->catfile($tempdir,'all.pseudomolecules.fasta');
my $all_gff3_tmp  = File::Spec->catfile($tempdir,'all.gff3');
my ($all_ctg_tgt,$all_pm_tgt,$all_gff3_tgt) = contig_file('all',1,$contigs_path);

my $all_ctg_out = Bio::SeqIO->new(-format => 'fasta', -file => ">$all_ctg_tmp");
my $all_pm_out  = Bio::SeqIO->new(-format => 'fasta', -file => ">$all_pm_tmp");

my @publish_ops;
foreach my $chrnum ( @chromosome_nums ) {

  eval {
    my $agp_file = agp_file($chrnum,1,$agp_path);
    my $agp_published_as = published_as($agp_file)
      or die "no agp file found for chromosome $chrnum in $agp_path\n";

    my ($contigs_tgt,$pm_tgt,$gff3_tgt) = contig_file($chrnum,1,$contigs_path);

    my @contigs = named_contigs( $chrnum, agp_file => $agp_published_as->{fullpath} );

    my $contigs_temp = File::Spec->catdir( $tempdir, "$chrnum.contigs");
    my $pm_temp = File::Spec->catdir( $tempdir, "$chrnum.pseudomolecule");
    my $gff3_temp = File::Spec->catdir( $tempdir, "$chrnum.gff3");

    # make a contigs fasta file
    { my $ctg_out = Bio::SeqIO->new( -format => 'fasta', -file => ">$contigs_temp");
      my @contigs_consume = @contigs;
      while ( my ($ctgname,$ctg) = splice @contigs_consume,0,2 ) {
	my $seq = agp_contig_seq($ctg, fetch_bac_sequence => \&fetch_bac_sequence, pad_short_sequences => 1)
	  or die "failed to build sequence for contig $ctgname";
	$seq = Bio::PrimarySeq->new( -id => $ctgname, -seq => $seq );
	$ctg_out->write_seq($seq);
	$all_ctg_out->write_seq($seq);
      }
    }

    # make a pseudomolecule fasta file
    my $pm_length = 0;
    { my $pm_out = Bio::SeqIO->new( -format => 'fasta', -file => ">$pm_temp");
      my $large_seq = agp_to_seq( $agp_published_as->{fullpath}, fetch_bac_sequence => \&fetch_bac_sequence, pad_short_sequences => 1 );
      $pm_length = $large_seq->length;
      my $bn = basename($agp_published_as->{fullpath});
      $large_seq->desc( "generated_from_agp_file:$bn" );
      $pm_out->write_seq( $large_seq );
      $all_pm_out->write_seq( $large_seq );
    }

    # make a gff3 file showing the composition of the chromosome
    { my @features = contig_features( @contigs );
      my $seqregion = Bio::SeqFeature::Annotated->new( -start => 1, -end => $pm_length, -seq_id => $features[0]->seq_id );
      my $gff3_out = Bio::FeatureIO->new( -format => 'gff', -version => 3, -file => ">$gff3_temp",
					  -sequence_region => $seqregion,
					);
      foreach ( @features ) {
	$gff3_out->write_feature( $_ );
      }
    }

    push @publish_ops,
      ['cp',$contigs_temp,$contigs_tgt],
      ['cp',$pm_temp,$pm_tgt],
      ['cp',$gff3_temp,$gff3_tgt];
  };
  if( $EVAL_ERROR ) {
    warn "could not build contigs for chromosome $chrnum:\n$EVAL_ERROR\n";
    next;
  }

}

$all_ctg_out = undef;
$all_pm_out  = undef;

# make the combined contigs gff3 file
CXGN::Tools::Run->run( 'gff3_reformat.pl',
		       -S => $all_pm_tmp,
		       (map $_->[1], grep $_->[2] =~ /\.gff3/, @publish_ops),
		       { out_file => $all_gff3_tmp },
		     );

push @publish_ops, ( ['cp',$all_ctg_tmp,$all_ctg_tgt],
		     ['cp',$all_pm_tmp,$all_pm_tgt],
		     ['cp',$all_gff3_tmp,$all_gff3_tgt],
		   );

#now gzip everything and add suffixes to the publish from and toe
foreach (@publish_ops) {
  my $r = CXGN::Tools::Run->run('gzip', '-c', '--rsyncable', '-n',  $_->[1]);
  CXGN::Tools::Run->run('cp', $r->out_file, $_->[1]);
}

#now do our publishing operations
publish(@publish_ops);

unlock_script();

#### SEQUENCE FETCHING ROUTINES ####
# each function is named type_fetch
# where type is the string returned by a call to CXGN::Tools::Identifiers::identifier_namespace

sub fetch_bac_sequence {
  my ($ident) = @_;
  return _get_chado_seq($ident);
}

#fetch a sequence by name from chado
sub _get_chado_seq {
  my ($ident) = @_;
  our $dbh ||= CXGN::DB::Connection->new;
  our $chado_feature_seq_q ||= $dbh->prepare('select residues from feature where name = ? limit 1');
  $chado_feature_seq_q->execute($ident)
    or return;
  my ($seq) = $chado_feature_seq_q->fetchrow_array;
  $chado_feature_seq_q->finish;
  return $seq;
}
