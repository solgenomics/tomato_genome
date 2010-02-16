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

    script_files => [ grep -f && -x, glob('scripts/*'), glob('scripts/old/*') ],

    # and now the list of perl module dependencies

    requires => {
        'autodie'            => 0,
        #'Bio::Root::Version' => '==1.005002',
        'Bio::Chado::Schema' => '0.05',

        'Cache::File' => 0,

        #'Chart::Lines' => 0,

        'DateTime'                   => 0,
        'DateTime::Format::Strptime' => 0,
        'DBI'                        => 0,
        'DBIx::Connector'            => 0,
        'Digest::MD5'                => 0,
        'enum'                       => 0,
        'File::Temp'                 => 0,
        'Hash::Util'                 => 0,
        'HTML::Entities'             => 0,
        'List::MoreUtils'            => 0,
        'List::Util'                 => 0,
        'LWP::Simple'                => 0,
        'Mail::Sendmail'             => 0,
        'Memoize'                    => 0,
        'Moose'                      => 0,
        'Moose::Role'                => 0,
        'MooseX::App::Cmd'           => 0,
        'MooseX::App::Cmd::Command'  => 0,
        'namespace::autoclean'       => 0,
        'Path::Class'                => 0,
        'Storable'                   => 0,
        'XML::DOM'                   => 0,
        'XML::LibXML'                => 0,
    },
    build_requires => {
        'aliased'           => 0,
        'App::Cmd::Tester'  => 0,
        'autodie'           => 0,
        'File::Temp'        => 0,
        'IO::String'        => 0,
        'IPC::Cmd'          => 0,
        'List::MoreUtils'   => 0,
        'List::Util'        => 0,
        'Path::Class'       => 0,
        'Test::Exception'   => 0,
        'Test::More'        => 0,
        'Test::Output'      => 0,
        'Test::Warn'        => 0,
        'Test::XML::Simple' => 0,
        'Tie::Function'     => 0,
    },

);
$build->create_build_script;
