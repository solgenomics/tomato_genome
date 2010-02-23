package CXGN::TomatoGenome::CmdLine::Command::load_gbrowse;
sub abstract {
    'update a versioned GBrowse DB from the given BAC repository'
}

use Moose;
use namespace::autoclean;

use CXGN::DB::GFF::Versioned;
use CXGN::TomatoGenome::BACPublish qw/ aggregate_filename publisher /;
use CXGN::TomatoGenome::Config;
use CXGN::IndexedLog;

extends 'CXGN::TomatoGenome::CmdLine::Command';
   with 'CXGN::TomatoGenome::CmdLine::DBConnector' => { connection_name => 'gbrowse_db' };
   with 'CXGN::TomatoGenome::CmdLine::DBConnector' => { connection_name => 'log_db'     };
   with 'CXGN::TomatoGenome::BACPublisher';

sub execute {
    my ( $self, $opt, $args ) = @_;

    my ($seqs,$gff3) =
        map { my $p = publisher()->published_as($_) or die "no version of $_ seems to be published!"; $p->{fullpath} }
        map aggregate_filename($_, $self->bac_repository_root),
        qw( all_seqs  all_gff3 );

    $self->gbrowse_db_conn;

    #have we loaded these already?
    my $log = CXGN::IndexedLog->open( DB => $self->log_db_conn->dbh, CXGN::TomatoGenome::Config->load_locked->{'bac_processing_log'} );
    my $log_string = "RELOADED_GBROWSE_WITH $gff3";


    unless( $log->lookup(content => $log_string) ) {

        my $bdb = CXGN::DB::GFF::Versioned->new( dsn      => $self->gbrowse_db_dsn,
                                                 user     => $self->gbrowse_db_user,
                                                 password => $self->gbrowse_db_password,
                                               );

        $self->vsay( "loading seqs: $seqs") if $seqs;
        $self->vsay( "loading gff3: $gff3") if $gff3;
        $bdb->load_new($seqs,$gff3);

        $log->append($log_string) unless $self->debug;
    }
}


__PACKAGE__->meta->make_immutable;
1;

