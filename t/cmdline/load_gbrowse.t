use strict;
use warnings;
use autodie ':all';

use FindBin;
use Test::More;
use File::Temp;
use File::Copy;
use Fatal qw/ copy /;

use Data::Dumper;

use aliased 'CXGN::TomatoGenome::CmdLine' => 'CmdLine';
use App::Cmd::Tester;
use Path::Class;

use CXGN::TomatoGenome::CmdLine;
use CXGN::TomatoGenome::CmdLine::Command::db_load;

use CXGN::DB::GFF::Versioned;

use lib './t/lib';
use FakeBACRepository;

my $test_dsn = $ENV{TOMATO_GENOME_TEST_DSN}
    or plan skip_all => 'must set TOMATO_GENOME_TEST_DSN to test load_gbrowse command', 2;

my $test_repos = FakeBACRepository->new;
my $gbrowse_test_db_root = 'test_load_gbrowse';
my $gbrowse_test_dsn = $test_dsn; $gbrowse_test_dsn =~ s/dbname=[^;]+/dbname=$gbrowse_test_db_root/;

# delete any existing test dbs
my $db = CXGN::DB::GFF::Versioned->new( dsn => $gbrowse_test_dsn );
$db->_clean_up_older_databases(0);

my $result = test_app( CmdLine, [ 'load_gbrowse',
				  '--debug',
				  '--bac_repository_root' => ''.$test_repos->bac_repository_root,
				  '--log_db_dsn'         => $test_dsn,
				  '--gbrowse_db_dsn' => $gbrowse_test_dsn,
				]);

like $result->stdout, qr(loading seqs:),
     'output mentions loading seqs'
    or diag Dumper $result;

like $result->stdout, qr(loading gff3:),
     'output mentions loading seqs'
    or diag Dumper $result;

is $result->error, undef, 'no error'
    or diag $result->error;

my $err = $result->stderr;
$err =~ s/^DBD::Pg[^\n]+\n//; #< ignore DBD::Pg warnings
is $err, '',
    'nothing on stderr'
    or diag Dumper $result;

$db = CXGN::DB::GFF::Versioned->new( dsn => $gbrowse_test_dsn );
like $db->_highest_db_version, qr/^\d+$/, 'seemed to have loaded something';

done_testing;
