package CXGN::TomatoGenome::CmdLine::DBConnector;
use MooseX::Role::Parameterized;
use DBIx::Connector;

parameter 'connection_name' => qw(
    is        ro
    isa       Str
    default   db
   );

parameter 'connection_description' => (
    is => 'ro',
    isa => 'Str',
    lazy_build => 1,
   ); __PACKAGE__->meta->parameters_metaclass->add_method(
   '_build_connection_description' => sub {
       my $n = shift->connection_name;
       $n =~ s/_/ /g;
       return $n;
   });


role {

    my $p    = shift;
    my $conn = $p->connection_name;
    my $desc = $p->connection_description;

    has $conn.'_dsn' => (
        documentation => "DBI dsn for connecting to the $desc",
        traits        => [qw(Getopt)],
        isa           => 'Str',
        is            => 'ro',
        required      => 1,
       );
    has $conn.'_user' => (
        documentation => "username for connecting to the $desc",
        traits        => [qw(Getopt)],
        isa           => 'Str',
        is            => 'ro',
       );
    has $conn.'_password' => (
        documentation => "password for connecting to the $desc",
        traits        => [qw(Getopt)],
        isa           => 'Str',
        is            => 'ro',
       );


    method $conn.'_conn' => sub {
        my ($self) = @_;

        no strict 'refs';
        $self->{$conn.'_dbc'} ||= DBIx::Connector->new( $self->{"${conn}_dsn"},
                                                        $self->{"${conn}_user"},
                                                        $self->{"${conn}_password"},
                                                      );
    };

};

# requires 'validate_args';
# requires 'usage_error';

# after    'validate_args' => sub {
#     my ( $self, $opt, $args ) = @_;

#     foreach my $required ('db_dsn','db_user','db_password') {
# 	$opt{$required} or $self->usage_error("must provide --$required option");
#     }
# };

1;
