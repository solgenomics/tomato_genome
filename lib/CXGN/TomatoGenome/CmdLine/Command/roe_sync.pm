package CXGN::TomatoGenome::CmdLine::Command::roe_sync;
sub abstract {
    'get Roe Lab BACs from GenBank and create .tar.gz BAC submissions for them' }

use Moose;
use namespace::autoclean;
use autodie qw/ :all /;

use Carp;

use LWP::Simple;
use Path::Class;
use File::Temp qw/ tempdir /;

use Bio::DB::GenBank;
use Bio::PrimarySeq;
use Bio::SeqIO;

use CXGN::TomatoGenome::BACPublish qw/ sequencing_files parse_filename /;

extends qw(CXGN::TomatoGenome::CmdLine::Command);

has 'table_url' => (
    documentation => 'the URL for the Roe table page',
    traits        => [qw(Getopt)],
    isa           => 'Str',
    is            => 'ro',
    default       => 'http://www.genome.ou.edu/tomato_table.html',
    cmd_aliases   => 'U',
);

has 'ftpsite_root' => (
    documentation => 'root directory of BACs ftp site',
    traits        => [qw(Getopt)],
    isa           => 'Str',
    is            => 'ro',
    required      => 1,
    cmd_aliases   => 'd',
);

has 'submission_destination' => (
    documentation => 'the destination path / URL to which new submissions will be copied',
    traits        => [qw(Getopt)],
    isa           => 'Str',
    is            => 'ro',
    default       => '/data/shared/tomato_genome/country_uploads/manual/upload',
    cmd_aliases   => 's',
);

has '_metadata' => (
    documentation => 'CXGN::Metadata object used for updating BAC metadata',
    is            => 'ro',
    isa           => 'CXGN::Metadata',
    lazy_build    => 1,
   ); sub _build__metadata { CXGN::Metadata->new }

has '_bac_status_log' => (
    documentation => 'CXGN::People::BACStatusLog object used for updating BAC metadata',
    is            => 'ro',
    isa           => 'CXGN::People::BACStatusLog',
    lazy_build    => 1,
   ); sub _build__bac_status_log {
       CXGN::People::BACStatusLog->new( shift->_metadata );
   }

sub execute {
    my ( $self, $opt, $args ) = @_;

    # get the roe html
    my $html = get( $self->table_url )
      or croak 'could not fetch tomato table from ' . $self->table_url;

    # parse it
    my $table = RoeTable->new( $html );

    my $gb_noseq = Bio::DB::GenBank->new( -retrievaltype => 'tempfile',
					  -seq_start => 1, -seq_stop  => 1,
					);
    my $gb_full  = Bio::DB::GenBank->new( -retrievaltype => 'tempfile' );

    # fetch the genbank record for it
    foreach my $table_rec (@{ $table->records }) {
	$self->vprint( "fetching $table_rec->{clone_name} / $table_rec->{gbacc}...\n" );
	my $gb_rec = RoeGenbankRecord->new( gb_richseq =>
						$gb_full->get_Seq_by_acc( $table_rec->{gbacc} )
						    || die( "no GenBank record found for $table_rec->{gbacc}" ),
					    table_rec => $table_rec,
					  );

	# if it has no clone, warn about it and skip
	unless( $gb_rec->clone ) {
	    warn "WARNING: could not find clone for accession ".$gb_rec->gb_richseq->accession.", skipping!\n";
	    next;
	}

	# if it does not have a chromosome assignment, give it one
        my $clone_reg_info = $gb_rec->clone->reg_info_hashref;
	unless( defined $clone_reg_info->{seq_proj}->{val} ) {
            unless( defined $gb_rec->chromosome_num ) {
                warn( $gb_rec->clone->clone_name." has no chromosome assignment in BAC registry, and cannot deduce it from its GenBank record.  Please manually set a chromosome number. Skipping.\n" );
                next;
            }

	    $self->vprint( 'Assigning '.$gb_rec->clone->clone_name.' to chromosome '.$gb_rec->chromosome_num."\n" );

            # look up the proper project id to use
            my ($proj_id) = $self->_metadata->selectrow_array( <<EOQ, undef, 'Tomato Chromosome '.$gb_rec->chromosome_num.' %' );
select sp_project_id from sp_project where name like ?
EOQ
            unless( $proj_id ) {
                warn "cannot find project ID for 'Tomato Chromosome ".$gb_rec->chromosome_num."'";
                next;
            }

            $self->_metadata->attribute_bac_to_project(
                $gb_rec->clone->clone_id,
                $proj_id
               );
            $self->_metadata->commit;

            # refresh the clone reg info
            $clone_reg_info = $gb_rec->clone->reg_info_hashref;
	}

        # if it has a chromosome assignment, but not a sequencing status, set it as 'complete'
        if( defined $clone_reg_info->{seq_proj}->{val} && $clone_reg_info->{seq_status}->{val} eq 'none' ) {
            my $current_proj = $self->_metadata->get_project_associated_with_bac( $gb_rec->clone->clone_id )
                or die "no current project??  this should not happen";

	    $self->vprint( 'Setting '.$gb_rec->clone->clone_name." sequencing status to 'complete'\n" );

            $self->_bac_status_log
                 ->change_status( bac        => $gb_rec->clone->clone_id,
                                  person     => 290, #< Robert Buels
                                  seq_status => 'complete',
                                 );
            $self->_bac_status_log->get_dbh->commit;
        }

	# check if we have the sequence on our ftp site, and that it is up to date
	my %files = sequencing_files( $gb_rec->clone, $self->ftpsite_root );
	next if $gb_rec->matches_seq_file( $files{seq} );
	
	# if not, make a BAC submission for it
	my $new_submission_file = $self->make_bac_submission( $gb_rec );

	# and copy it into place for submission to the pipeline
	-d $self->submission_destination or croak "dir ".$self->submission_destination." does not exist!";
	my @cmd = ( 'cp', $new_submission_file, $self->submission_destination );
	system @cmd;
    }
}

# makes a submission tarball for the given object and returns a Path::Class::File
sub make_bac_submission {
    my ( $self, $gb_rec ) = @_;

    # make a BAC submission tarball for it
    my $clone_name_with_chrom = $gb_rec->clone->clone_name_with_chromosome
	or die $gb_rec->clone->clone_name." has no chromosome assignment!\n";

    my $tempdir = tempdir( CLEANUP => 1 ); #< will last until program exit
    my $bac_dir = dir( $tempdir, $clone_name_with_chrom );

    $bac_dir->mkpath or die "could not make dir '$bac_dir'";

    # write the genbank accession
    $bac_dir->file('gbacc.txt')->openw->print($gb_rec->accession."\n");

    my $seq = Bio::PrimarySeq->new( -id   => $clone_name_with_chrom,
				    -seq  => $gb_rec->seq,
				   );
    # write the sequence file
    Bio::SeqIO->new( -format => 'fasta',
		     -file => '>'.$bac_dir->file("$clone_name_with_chrom.seq"),
		    )
	      ->write_seq( $seq );
    my $tarball = file($tempdir, "$clone_name_with_chrom.tar.gz");
    utime 0, 0, $bac_dir, $bac_dir->children; #< set the file
    system "tar -C $tempdir -m -cf - $clone_name_with_chrom | gzip -c -n --rsyncable > $tarball";

    $bac_dir->rmtree; #< clean up some of the temp files

    return $tarball;
}

# class representing a Roe-style HTML table of stuff
package RoeTable;
use Moose;

has 'records' => (
    isa => 'ArrayRef',
    is  => 'ro',
    required => 1,
);

around BUILDARGS => sub {
    my ( $orig, $class, $html ) = @_;
    my $recs = $class->parse_roe_table( $html );
    return $class->$orig( records => $recs );
};

# takes roe html, returns arrayref of sequence records
sub parse_roe_table {
    my ( $class, $html ) = @_;
    die 'must pass html' unless $html;

    my @recs;
    open my $f, '<', \$html;
    while ( my $line = <$f> ) {
        chomp $line;
        $line =~ s/^\s*|\s*$//g;

        no warnings 'uninitialized';

        my %p;
        ( $p{type}, $p{clone_name} ) = split /\s+/, $line;
        ( $p{gbacc} ) = $line =~ /([A-Z]{2}\d{6,7})/
          or next;
        ( $p{htgs_phase} ) = $line =~ /$p{gbacc}.{0,7}\s+(\d+)/
          or next;

        push @recs, \%p;
    }

    return [ reverse @recs ];
    return \@recs;
}


__PACKAGE__->meta->make_immutable;

# class representing a Roe-style genbank record, encapsulating a
# Bio::Seq::RichSeq
package RoeGenbankRecord;
use Moose;
use CXGN::Genomic::CloneIdentifiers ();

has 'gb_richseq' => (
    isa => 'Bio::Seq::RichSeq',
    is => 'rw',
    required => 1,
    handles => [ 'seq', 'length', 'accession' ],
);
has 'table_rec'  => (
    isa => 'HashRef',
    is => 'ro',
    required => 1,
);

has 'clone' => (
     is => 'ro',
     lazy_build => 1,
);
sub _build_clone {
    my $self = shift;
    my $seq = $self->gb_richseq;
    my $clone_name = $self->clone_name
        or return;

    require CXGN::Genomic::Clone;
    return CXGN::Genomic::Clone->retrieve_from_clone_name( $clone_name );
}


has 'clone_name' => (
    is => 'ro',
    lazy_build => 1,
);
sub _build_clone_name {
    my $self = shift;
    my $seq = $self->gb_richseq;

    my @tests =
        (
         sub { shift->desc =~ /clone\s+([\w\-]+)/ && $1 },
         sub { my $a = shift->annotation(); my ($c) = $a->get_Annotations('clone'); $c && $c->value },
        );

    foreach my $t (@tests) {
        my $name = $t->($seq);
        return $name if $name && CXGN::Genomic::CloneIdentifiers::parse_clone_ident($name);
    }

    return;
}


has 'chromosome_num' => (
    is => 'ro',
    lazy_build => 1,
   ); sub _build_chromosome_num {
       shift->gb_richseq->desc =~ /chromosome (\d+)/i
           and return $1;
       return;
   }

# check whether this genbank record matches the sequence in the given file
sub matches_seq_file {
    my ( $self, $seqfile ) = @_;

    return unless $seqfile && -f $seqfile;

    my $seq = Bio::SeqIO->new( -format => 'fasta', -file => $seqfile )->next_seq
	or return;

    return unless $seq->length == $self->length
	&& $seq->seq eq $self->seq;
    return 1;
}

__PACKAGE__->meta->make_immutable;

###
1;#
###
