use strict;
use warnings;
use Carp;

use Test::More tests => 10;
use Test::Output qw/ output_from /;
use App::Cmd::Tester;

use FindBin;
use Path::Class;

use Data::Dumper;

use File::Temp;

use Bio::DB::GenBank;
use Bio::SeqIO;

use_ok('CXGN::TomatoGenome::CmdLine');
use_ok('CXGN::TomatoGenome::CmdLine::Command::roe_sync');

# we have a test copy of the tomato table here
my $test_tomato_html_file = dir( $FindBin::RealBin )->parent->subdir('data')->file( 'tomato_table.html' );

my $table = RoeTable->new( scalar $test_tomato_html_file->slurp );
my $test_rec = $table->records->[0];
is_deeply( $test_rec,
	   {
	       'htgs_phase' => '3',
	       'gbacc' => 'AC238499',
	       'type' => 'BAC',
	       'clone_name' => 'hba-8l19'
	   },
	   'parsed OK' )
    or diag Dumper($table->records->[0]);

my $gb_full  = Bio::DB::GenBank->new( -retrievaltype => 'tempfile' );

# get a seq from genbank
my $test_seq = $gb_full->get_Seq_by_acc( $test_rec->{gbacc} );

# write it to a tempfile
my $test_seq_file = File::Temp->new;
Bio::SeqIO->new( -fh => $test_seq_file, -format => 'fasta' )->write_seq( $test_seq );

my $gbrec = RoeGenbankRecord->new( gb_richseq => $test_seq,
				   table_rec  => $test_rec,
				  );

# test the matches_seq_file method
ok( $gbrec->matches_seq_file( "$test_seq_file" ), 'matches_seq_file method works 1' );
$gbrec->gb_richseq->seq('ACTCGTACGA');
ok( ! $gbrec->matches_seq_file( "$test_seq_file" ), 'matches_seq_file method works 2' );


# test make_bac_submission
my $file = CXGN::TomatoGenome::CmdLine::Command::roe_sync::make_bac_submission( undef, $gbrec );
isa_ok( $file, 'Path::Class::File' );
ok( -f $file, 'bac submission file exists')
    or diag $file;
my $list = `tar -tzf $file`;
is( $list, <<EOF, 'new submission tarball contains the right stuff' );
C01HBa0008L19/
C01HBa0008L19/C01HBa0008L19.seq
C01HBa0008L19/gbacc.txt
EOF

# now run integration tests
SKIP: {
    skip 'set CXGN_ROE_SYNC_LIVE_GENBANK=1 to run live genbank tests', 2
	unless $ENV{CXGN_ROE_SYNC_LIVE_GENBANK};

    my $cmd = CXGN::TomatoGenome::CmdLine->new;

    my $tempdir = File::Temp->newdir;
    my $destination = dir( $tempdir, 'test_destination' );
    $destination->mkpath;

    local @ARGV = ( 'roe_sync',
		    '--table_url'    => 'file://'.$test_tomato_html_file,
		      '--ftpsite_root' => "$tempdir",
		    '--submission_destination' => "$destination",
		   );
    my ($output,$err) = output_from( sub{ $cmd->run } );
    like $output,
	 qr/fetching/,
	 'output looks OK';

    ok( -f $destination->file( $file->basename ), 'submission file was put in place' )
	  or diag $output.`find $destination`;
}
