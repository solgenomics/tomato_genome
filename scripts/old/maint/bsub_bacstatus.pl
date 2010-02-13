#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use File::Basename;
use File::Spec;

use CXGN::TomatoGenome::BACSubmission;
use CXGN::Genomic::CloneIdentifiers qw/ parse_clone_ident /;
use CXGN::Genomic::Clone;

my $conf = CXGN::TomatoGenome::Config->load_locked;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script <bac name> <bac name> ...

  Given a list of BACs on the command line, check up on each of them

  Options:

    -V skip validation of any newly uploaded tarballs

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('V',\%opt) or usage();

my @reports;
foreach my $bacname (@ARGV) {
  @reports = ();

  #check that the clone name is parsable
  my $p = parse_clone_ident($bacname);
  unless( $p ) {
    push @reports, "ERROR: invalid clone identifier, please check its format";
    next;
  }

  #check that it's in the database
  my $clone = CXGN::Genomic::Clone->retrieve_from_parsed_name($p);
  unless( $clone ) {
    push @reports, "ERROR: clone not found in database";
  }

  if ( $clone ) {
    #check for any uploaded files
    if ( my $clone_name = $clone->clone_name_with_chromosome ) {
      my $globexpr = File::Spec->catfile($conf->{'country_uploads_path'},'*','upload',"$clone_name*.tar.gz");
      my @uploaded_files = glob $globexpr;
      if ( @uploaded_files ) {
	foreach my $ulfile (@uploaded_files) {
	  
	  #parse out the account that uploaded it
	  my @dirs = File::Spec->splitdir($ulfile);
	  my $bn = pop @dirs;
	  pop @dirs;		#< upload
	  my $account = pop @dirs;

	  my $ul_report =  "uploaded file still in pipeline - account:$account  file:upload/$bn";

	  unless ( $opt{V} ) {
	    eval {
	      my $sub = CXGN::TomatoGenome::BACSubmission->open($ulfile);

	      if (my @validation_errors = $sub->validation_errors) {
		$ul_report .= "\n      VALIDATION ERRORS:";
		$ul_report .= "\n        - ".$sub->error_string($_) foreach @validation_errors;
	      } else {
		$ul_report .= " passes validation";
	      }
	    }; if($EVAL_ERROR) {
	      push @reports, "ERROR: $EVAL_ERROR\nCould not open uploaded file, perhaps it is unfinished or corrupt?";
	    }
	  }
	  push @reports, $ul_report;
	}
      }
    } else {
      push @reports, 'no uploaded tarfiles stuck in pipeline';
    }
    
  } else {
    push @reports, "ERROR: not assigned to sequencing project in SGN BAC registry";
  }
  
  #check its BAC registry status
  my $reg = $clone->reg_info_hashref;
  push @reports, "SGN BAC registry contains: (seq chr: $reg->{seq_proj}->{disp}, seq status: $reg->{seq_status}->{disp})";
  unless( $reg->{seq_proj}->{disp} eq $p->{chr} ) {
    my $chrstr = 'none' if $reg->{seq_proj}->{disp} eq '-';
    push @reports, "ERROR: SGN BAC registry has chromosome '$chrstr', but BAC name indicates chromosome '$p->{chr}'";
  }
  
  #check for any published files
  my %files = ftp_sequencing_files($clone);
  if ($files{tar}) {
    push @reports, "published tarfile $files{tar}";
  } else {
    push @reports, "no tarfile passed";
  }
  if ($files{seq}) {
    push @reports, "published seq file $files{seq}";
  } else {
    push @reports, "no sequence passed";
  }

  #check for loaded sequence
  if ( my @seqnames = $clone->latest_sequence_name ) {
    if ( @seqnames > 1 ) {
      push @reports, "ERROR: sequence is in multiple fragments, please re-upload to conform with SOL Bioinformatics Guidelines";
    }

    push @reports, "sequences found in SGN DB: (".join(',',@seqnames).")";
    my $sp = $clone->seqprops;
    if ( $sp && $sp->{htgs_phase} ) {
      push @reports, "HTGS phase in SGN DB: $sp->{htgs_phase}";
    } else {
      push @reports, "ERROR: no HTGS phase recorded in SGN DB";
    }

    #check if it has a genbank accession
    if (my $acc =  $clone->genbank_accession ) {
      push @reports, "recorded GenBank accession in SGN DB: $acc";
    } else {
      push @reports, "ERROR: no valid GenBank accession recorded in SGN DB"
    }
  } else {
    push @reports, "no sequences in SGN DB\n";
  }

} continue {
  print "\n== $bacname ===========================\n";
  foreach (@reports) {
    print "   - $_\n";
  }
}


=head2 ftp_sequencing_files

  Usage: my %files = ftp_sequencing_files($clone)
  Desc : get ftp links for annotation and other files related to this clone
  Ret  : a hash-style list as:
         (  seq  => 'ftp://ftp.sgn.cornell.edu/tomato_genome/bacs/chrXX/CXXblahblah.seq',
            tar  => the tarfile,
            gff3    => the gff3 file of all automatic annotations to this clone,
            gamexml => the gamexml file of all automatic annotations to this clone,
         )
         or undef if the files can't be found
  Args : none
  Side Effects: looks things up in the filesystem

=cut

sub ftp_sequencing_files {
  my ($self) = @_;

  my $ftp_path = $conf->{'ftpsite_root'}
    or die "ftpsite_root configuration variable is not defined";
  my $ftp_url = $conf->{'ftpsite_url'}
    or die "ftpsite_url configuration variable is not defined";

  my %files = $self->sequencing_files($ftp_path);

  #replace the ftpsite's filesystem root with the ftp site url
  foreach(values %files) {
    $_ =~ s/^$ftp_path/$ftp_url/e if $_;
  }

  return %files;
}
