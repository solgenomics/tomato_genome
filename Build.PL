use strict;
use Module::Build;

my $build = Module::Build->new(
    dist_name     => 'Tomato-Genome',
    dist_author   => 'Sol Genomics Network',
    dist_abstract => 'scripts and libraries for managing tomato genome data',
    license       => 'perl',
    create_makefile_pl => 'passthrough',

    # current version of our core libraries distribution
    dist_version => '0.1',

    recursive_test_files => 1,

    # and now the list of perl module dependencies
    requires => {

        # right now, we have to depend on whole bioperl, like this
        'Bio::Root::Version' => '1.005002',

        #'Bio::DB::GenBank'                => 0,

        #'CXGN::Genomic::Clone'            => 0,
        #'CXGN::Genomic::CloneIdentifiers' => 0,

	'Bio::Chado::Schema'           => '0.05',
	'DBIx::Connector'              => 0,

        'LWP::Simple'                  => 0,
        'Moose'                        => '0.93',
        'MooseX::App::Cmd'             => '0.06',
        'namespace::autoclean'         => '0.09',
    },
    build_requires => {
	'aliased'          => 0,
        'App::Cmd::Tester' => 0,
        'Data::Dumper'     => 0,
        'File::Temp'       => 0,
        'Fatal'            => 0,
        'FindBin'          => 0,
        'Path::Class'      => 0,
	'Test::Output'     => 0,
    },
);
$build->create_build_script;

