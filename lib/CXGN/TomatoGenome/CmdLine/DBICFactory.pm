package CXGN::TomatoGenome::CmdLine::DBICFactory;
use Moose::Role;
use Carp ();

requires 'db_dsn', 'db_user', 'db_password';

sub open_dbic_schema {
    my ($self, $package_name, %options) = @_;

    eval "require $package_name";
    Carp::confess "could not require $package_name: $@" if $@;

    # coerce search path to appropriate string
    $options{search_path} = join(',',@{$options{search_path}}) if ref $options{search_path};

    return $package_name
             ->connect(  $self->db_dsn,
			 $self->db_user,
			 $self->db_password,
			 $options{dbconn_args} || {},
			 {
                          ( $options{search_path}
                              ? ( on_connect_do => [ 'SET search_path TO '.$options{search_path} ] )
                              : ()
			  )
			 },
		       );

}

###
1;#
###

