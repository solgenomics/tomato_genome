use strict;
use warnings;
use Test::More;
use Test::Warn;

use FindBin;
use Path::Class;

use_ok('CXGN::DB::GFF::Versioned');

# test compat mode
{ my $d;
  warning_like {
      $d = CXGN::DB::GFF::Versioned->new( -db => 'db_gff_versioned_test' );
  } qr/deprecated/, 'got deprecation warning';

  is $d->dbname(231), 'db_gff_versioned_test.231', 'good dbname';
}


SKIP: {
    my $test_dsn = $ENV{TOMATO_GENOME_TEST_DSN};
    skip 'set TOMATO_GENOME_TEST_DSN to enable database tests'
        unless $test_dsn;

    $test_dsn =~ s/dbname=[^;]+/dbname=db_gff_versioned_test/
        or die "no dbname found in test dsn '$test_dsn'";

    my $db = CXGN::DB::GFF::Versioned->new( dsn => $test_dsn );

    $db->_clean_up_older_databases(0);

    is $db->_highest_db_version, undef, 'correct highest db version';

    my $seq_file = file( $FindBin::RealBin, 'data', 'C06HBa0012B01.1.v87.seq' );
    $db->load_new([$seq_file]);
    is $db->_highest_db_version, 1, 'correct highest db version';
    $db->load_new([$seq_file]);
    is $db->_highest_db_version, 2, 'correct highest db version';

    $db->_clean_up_older_databases(0);

    is $db->_highest_db_version, undef, 'correct highest db version';
}

done_testing;

###########################

sub opendb {
    CXGN::DB::GFF::Versioned->new( dsn => $ENV{TOMATO_GENOME_TEST_DSN} );
}
