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

$ENV{TOMATO_GENOME_TEST_DSN}
    or plan skip_all => 'must set TOMATO_GENOME_TEST_DSN to test db_load command', 2;


my $tempdir = File::Temp->newdir;
my $tf = 'C06HBa0012B01.1.v87.seq';
my $test_file = dir( $FindBin::RealBin )->parent->subdir('data')->file($tf);
my $dest_file = dir($tempdir)->subdir('chr06')->subdir('finished')->file($tf);
$dest_file->dir->mkpath;
copy( "$test_file", "$dest_file" );

my $result = test_app( CmdLine, [ 'db_load',
				  '--debug',
				  '--ftpsite_root' => "$tempdir",
				  '--db_dsn' => $ENV{TOMATO_GENOME_TEST_DSN},
				]);

like $result->stdout, qr(loading seq file .*$tf),
     'output mentions loading test file'
    or diag Dumper $result;

my $err = $result->stderr;
$err =~ s/^DBD::Pg[^\n]+\n//; #< ignore DBD::Pg warnings
is   $err, '',
     'nothing on stderr'
    or diag Dumper $result;

done_testing;
