package CXGN::TomatoGenome::CmdLine::Command::db_load;
use Moose;
use namespace::autoclean;
use Carp;

extends 'CXGN::TomatoGenome::CmdLine::Command';
   with 'CXGN::TomatoGenome::CmdLine::DBConnector';
   with 'CXGN::TomatoGenome::CmdLine::DBICFactory';

use Path::Class;
use Hash::Util qw/ lock_keys /;

use CXGN::Genomic::Clone;
use CXGN::TomatoGenome::BACPublish qw/ glob_pattern parse_filename /;

has 'ftpsite_root' => (
    documentation => 'root directory of BACs ftp site',
    traits        => [qw(Getopt)],
    isa           => 'Str',
    is            => 'ro',
    required      => 1,
    cmd_aliases   => 'F',
);

sub chado {
    my ($self) = @_;
    return $self->{chado} ||= $self->open_dbic_schema('Bio::Chado::Schema', search_path => 'public,genomic');
}

sub abstract {
    'load BACs from the file-based BAC repository into the DB'
}

sub validate_args {
    my ($self, $opt, $args) = @_;
    @$args and $self->usage_error("unknown argument '$args->[0]'");
}

sub execute {
    my ( $self, $opt, $args ) = @_;

    my $glob_pat = glob_pattern('all_seqs',$self->ftpsite_root );
    $self->vprint("searching for seqfiles matching $glob_pat...\n");
    foreach my $bac_seq ( glob glob_pattern('all_seqs',$self->ftpsite_root ) ) {
	if( $self->already_loaded( $bac_seq ) ) {
	    $self->vprint("already loaded $bac_seq\n");
	}
	else {
	    $self->vprint("loading seq file $bac_seq\n");
	    $self->load_bac_seq( $bac_seq );
	}
    }
}

sub already_loaded {
    my ( $self, $seqfile ) = @_;

    my $p = parse_filename( $seqfile )
	or confess "could not parse filename '$seqfile'";

    # check whether a feature for this seq exists
    return $self->tomato_organism
	        ->search_related('features', { name => $p->{seq_name} } )
		->count;
}

sub load_bac_seq {
    my ( $self, $seqfile ) = @_;

    my $seq = Bio::SeqIO->new( -file => $seqfile, -format => 'fasta' )->next_seq;

    my $desc = $seq->desc;
    my $dl_info = $self->parse_desc( $desc );

    my $clone = CXGN::Genomic::Clone->retrieve_from_clone_name( $seq->id )
	or die "could not retrieve clone for seq name ".$seq->id;

    # load it into the DB with the proper accession
    # and update the public.clone_feature table
    $self->chado->txn_do( sub {

        my $gbacc_version = delete $dl_info->{gb_version};
        my $gbacc         = delete $dl_info->{gb_accession};

        # make a feature for it in the feature table
        my $new_feature =
	      $self->tomato_organism
                   ->create_related('features',
                                    { name => $seq->id,
                                      uniquename => $seq->id,
                                      type  => $self->bac_cvterm,
                                      residues => $seq->seq,
                                      seqlen => $seq->length,
                                    },
                                   );
        $new_feature->create_featureprops({ %$dl_info,
                                            description => $desc,
                                            chromosome => $clone->chromosome_num || 0,
                                          },
                                          { autocreate => 1 },
                                         );

        # add a dbxref for its genbank accession
	if( $gbacc ) {
	    my $gbacc_dbx =
		$self->chado
		     ->resultset('General::Db')
		     ->find_or_create({ name => 'DB:GenBank_Accession'})
		     ->find_or_create_related('dbxrefs',
					      { accession => $gbacc,
						version   => $gbacc_version,
					      }
					     );
	    $new_feature->add_to_secondary_dbxrefs( $gbacc_dbx );
	}

        # manually update the clone_feature table
        $self->chado->storage->dbh_do(sub {
					  my ($s,$dbh) = @_;
					  $dbh->do('delete from clone_feature where clone_id = ?',
						   undef,
						   $clone->clone_id,
						  );
					  $dbh->do('insert into clone_feature (feature_id,clone_id) values (?,?)',
						   undef,
						   $new_feature->feature_id,
						   $clone->clone_id,
						  );
				      });

        # don't actually do anything to the db if debug is on
        $self->chado->storage->txn_rollback if $self->debug;
    });
}

# given a desc line containing some metadata, return a hashref with
# info parsed out of it, with keys gb_accession, gb_version,
# htgs_phase, and sequenced_by, not all of which are necessarily
# present
sub parse_desc {
  my ( $self, $desc ) = @_;

  lock_keys( my %info, qw/ htgs_phase gb_accession gb_version sequenced_by upload_account_name /);

  if($desc =~ /htgs_phase:(\d)/) {
    $info{htgs_phase} = $1;
    die "invalid phase $info{htgs_phase}" unless grep $info{htgs_phase} == $_, 1..3
  }

  if($desc =~ /^[A-Z_]{2,3}\d+\.(\d+)/) {
    $info{gb_accession} = $&;
    $info{gb_version}   = $1;
  }

  if( $desc =~ /sequenced_by:(\S+)/ ) {
    $info{sequenced_by} = $1;
  }

  if( $desc =~ /upload_account_name:(\S+)/ ) {
    $info{upload_account_name} = $1;
  }

  return \%info;
}

sub tomato_organism {
    shift
	->chado
        ->resultset('Organism::Organism')
        ->find_or_create({ genus => 'Solanum', species => 'lycopersicum' });
}

# get the cvterm for BAC clones from chado
sub bac_cvterm {
    shift
	->chado
        ->resultset('Cv::Cv')
        ->search({'me.name' => 'sequence'})
        ->search_related(cvterms => {'cvterms.name' => 'BAC_clone'})
        ->first;
}


__PACKAGE__->meta->make_immutable;
1;
