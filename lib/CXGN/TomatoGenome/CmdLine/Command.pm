package CXGN::TomatoGenome::CmdLine::Command;
use Moose;
extends qw/ MooseX::App::Cmd::Command /;

has 'quiet' => (
    documentation => 'turn off all output except errors',
    traits        => [qw(Getopt)],
    isa           => 'Bool',
    is            => 'ro',
    cmd_aliases   => 'q',
    default       => 0,
);

has 'debug' => (
    documentation => 'debug/test mode, developers only',
    traits        => [qw(Getopt)],
    isa           => 'Bool',
    is            => 'ro',
    default       => 0,
);

sub vprint {
    my $self = shift;
    print @_ unless $self->quiet && !$self->debug;
}

sub vsay {
    shift->vprint(@_,"\n");
}

__PACKAGE__->meta->make_immutable;

1;
