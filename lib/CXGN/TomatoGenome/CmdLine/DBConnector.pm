package CXGN::TomatoGenome::CmdLine::DBConnector;
use Moose::Role;
use DBIx::Connector;

has 'db_dsn' => (
    documentation => 'DBI Data Source Name (DSN) for the database to connect to',
    traits        => [qw(Getopt)],
    isa           => 'Str',
    is            => 'ro',
    required      => 1,
    cmd_aliases   => 'd',
);
has 'db_user' => (
    documentation => 'root directory of BACs ftp site',
    traits        => [qw(Getopt)],
    isa           => 'Str',
    is            => 'ro',
    cmd_aliases   => 'u',
);
has 'db_password' => (
    documentation => 'root directory of BACs ftp site',
    traits        => [qw(Getopt)],
    isa           => 'Str',
    is            => 'ro',
    cmd_aliases   => 'p',
);

sub dbc {
    my ($self) = @_;
    $self->{dbc} ||= DBIx::Connector->new( $self->db_dsn, $self->db_user, $self->db_password );
}

# requires 'validate_args';
# requires 'usage_error';

# after    'validate_args' => sub {
#     my ( $self, $opt, $args ) = @_;

#     foreach my $required ('db_dsn','db_user','db_password') {
# 	$opt{$required} or $self->usage_error("must provide --$required option");
#     }
# };

1;
