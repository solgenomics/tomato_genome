package Bio::DB::GFF::Adaptor::dbi::cxgn_pg;

use strict;
use warnings;
use base qw/Bio::DB::GFF::Adaptor::dbi::pg/;
our @ISA = ('Bio::DB::GFF::Adaptor::dbi::pg');

use Bio::DB::GFF::Util::Rearrange; # for rearrange()


use CXGN::DB::Connection;

use CXGN::DB::GFF::Versioned;

use CXGN::CDBI::Class::DBI;

use CXGN::VHost;


=head1 OVERVIEW

This is a subclass of L<Bio::DB::GFF::Adaptor::dbi::pg>, with modifications to make it use CXGN::DB::GFF::Versioned for connecting to the most recent version of the GFF database.

=cut

sub new {
  my $class = shift;
  my ($dbname_root) = rearrange([
				[qw(FEATUREDB DB DSN)],
				],@_);

  $dbname_root
    or die "must pass -db option to cxgn_pg database adaptor!\n";


  my $db = CXGN::DB::GFF::Versioned->new( -db => $dbname_root );
  my $bdb = $db->bdb
    or die "no current version of '$dbname_root' DB exists!";
  return $bdb;
}


###
1;#do not remove
###
