package CXGN::TomatoGenome::BACPublisher;
use Moose::Role;
use namespace::autoclean;

has 'bac_repository_root' => (
    documentation => 'root directory of BACs file repository',
    isa           => 'Str',
    is            => 'ro',
    required      => 1,
);


1;

