package CXGN::DB::GFF::Versioned;
use Moose;
use namespace::autoclean;
use Carp;

use Bio::DB::GFF;

=head1 NAME

CXGN::DB::GFF::Versioned - wrapper for Bio::DB::GFF with a Postgres
backend, providing versioned loading and cleanup of Bio::DB::GFF
databases

=head1 SYNOPSIS

  # open a handle on the versioned series of databases
  # named my_db_base_name.1, my_db_base_name.2, ...
  my $db = CXGN::DB::GFF::Versioned->new( dsn => $dsn, user => $user, password => $password );

  # load a new version of this DB with the given seq and gff3 file
  # will make my_db_base_name.(n+1) where n is the largest version of
  # my_db_base_name currently present.
  $db->load_new( $seq_file, $gff3_file );

  # get a Bio::DB::GFF object opened to this most recent db
  my $bdb = $db->bdb;

  # get a Bio::DB::GFF object opened to version 2 of this db.  if
  # there is no version 2, $bdb will be undef
  my $bdb = $db->bdb(2);


=head1 ATTRIBUTES

dsn - database DSN, required

user - database user, required

password - database password, required

=cut

for (qw( dsn user password )) {
    has $_ => (
        is  => 'ro',
        isa => 'Str',
       );
}
has '+dsn' => (required => 1);

=head1 METHODS

=head2 new

  Usage: my $dbgff = CXGN::DB::GFF::Versioned->new(
               dsn => $dsn,
               user => $user,
               password => $password
         );
  Desc : opens the most recent version of the given database name
  Args : hash-style list as:
           dsn      => DBI connection string for the database,
           user     => (optional) username to connect as, defaults to blank,
           password => (optional) password to use for connecting, defaults to null
  Ret  : new object

=cut

around BUILDARGS => sub {
    my $orig = shift;
    my $class = shift;
    my %args = @_;

    # did we get a -db arg?
    my $compat_dbname = delete $args{-db};

    # if it's the new interface, this is all we do
    return $class->$orig(@_) unless $compat_dbname;

    #### otherwise, if we got -db, have to do some things for compat
    warn "WARNING: -db option to CXGN::DB::GFF::Versioned is deprecated, pass dsn, user, and password instead";

    die "cannot specify both -db and other args for CXGN::DB::GFF::Versioned"
        if %args;

    require CXGN::DB::Connection;
    my ($dsn, $user, $pass) = CXGN::DB::Connection->new_no_connect
                                                  ->get_connection_parameters;

    if( $compat_dbname =~ /^cxgn_conf:(.+)$/ ) {
        require SGN::Config;
        $compat_dbname = SGN::Config->load_locked->{$1};
    }

    $dsn = $class->_replace_dbname( $dsn, $compat_dbname );

    return $class->$orig( dsn => $dsn, user => $user, password => $pass );

};


=head2 dbname

  Usage: my $n = $obj->dbname;
  Desc : get the (versioned) name of the database being used,
         e.g. 'bio_db_gff.2' for version 2 of the bio_db_gff
         database
  Args : optional version number, defaults to the latest
         version
  Ret  : text string, or nothing if no version number passed and
         no current database exists
  Side Effects: none

  If you pass an explicit version number, does not check whether that
  database actually exists.

=cut

sub dbname {
  my ($self,$version) = @_;
  $version ||= $self->_highest_db_version
    or return;
  $version =~ /\D/ && $version ne 'tmp' and confess "invalid version '$version'";
  $self->dbname_root or confess "no dbname_root defined";
  return $self->dbname_root.".$version";
}

has 'dbname_root' => ( is => 'ro',
                       lazy_build => 1,
                      ); sub _build_dbname_root {
                          my $d = shift->dsn;
                          my ($n) = $d =~ /dbname=([^;]+)/
                              or die "could not parse dbname from dsn $d";
                          $n =~ s/\.\d+$//;
                          return $n;
                      };

# dbh to use for creating databases, etc.  defaults to highest
# existing version of the db, or postgres if there is none
has 'control_dbh' => ( is => 'ro',
                       lazy_build => 1,
                      ); sub _build_control_dbh {
                          my $self = shift;
                          return DBI->connect( $self->_replace_dbname( $self->dsn, 'postgres' ),
                                               $self->user,
                                               $self->password,
                                              );
                      }

=head2 bdb

  Usage: my $bdb = $obj->bdb
  Desc : get a Bio::DB::GFF handle for one of the versions of the database
  Args : optional numerical version to open, defaults to most recent
  Ret  : Bio::DB::GFF object at the requested version
         or nothing if it does not exist
  Side Effects: may open a new database connection
  Example:

=cut

sub bdb {
  my ($self,$version) = @_;
  $version = $self->_highest_db_version() unless defined $version;
  return unless $version;
  #warn "bdb with version $version\n";

  return $self->_open_bdb($version);
}

sub _versioned_dsn {
    my ($self, $version) = @_;
    my $dbname = $self->dbname($version);
    return $self->_replace_dbname( $self->dsn, $dbname );
}
sub _replace_dbname {
    my ( $self, $dsn, $dbname ) = @_;
    $dsn =~ s/dbname=[^;]+/dbname=$dbname/i
        or confess "cannot parse dbname from dsn '$dsn'";
    return $dsn;
}

#open a bdb with the given database name, and possibly make it writable
sub _open_bdb {
  my ($self,$version,$new_flag) = @_;

  my $dbname = $self->dbname($version);

  my $dsn = $self->_versioned_dsn( $version );

  my $bdb= Bio::DB::GFF->new( -adaptor => 'dbi::pg',
			      -dsn     => $dsn,
                              $new_flag ? (-write => 1,-create => 1) : (),
                              -user    => $self->user,
                              -pass    => $self->password,
			    )
    or die "Can't open Bio::DB::GFF $dbname database: ",Bio::DB::GFF->error,"\n";

  return $bdb;
}


=head2 load_new

  Usage: $obj->load_new([$seq_file,$seq_file2,...],[$gff3_file, $gff3_file2, ...]);
  Desc : loads the given sequence and gff3 file into a new version of
         the database, and updates this handle to point to the new version
  Args : base name,
         arrayref of seq files or single seq file,
         arrayref of GFF3 files or single GFF3 file
  Ret  : new handle
  Side Effects: creates new databases, loads data into them, dies on
                load errors

=cut

sub load_new {
  my ($self,$seqs,$gff3) = @_;

  $seqs ||= [];
  $gff3 ||= [];
  $seqs = [$seqs] unless ref $seqs;
  $gff3 = [$gff3] unless ref $gff3;

  # make the .tmp database if necessary
  my $tmp_db = $self->dbname('tmp');
  $self->_make_db($tmp_db) unless $self->_db_exists($tmp_db);

  #open the gff and fasta files at the same time, so they don't get moved out from under us
  my @gff3_fh  = map { open my $f, $_ or die("$! opening '$_' for reading\n"); $f } @$gff3;
  my @fasta_fh = map { open my $f, $_ or die("$! opening '$_' for reading\n"); $f } @$seqs;

  # open our filehandles to /dev/null to shut up the idiotic warnings
  # and status messages spewed by this bioperl code
  local *STDOUT_SAVE;
  local *STDERR_SAVE;
  open STDOUT_SAVE, ">&STDOUT" or die "$! saving STDOUT";
  open STDERR_SAVE, ">&STDERR" or die "$! saving STDERR";
  open STDOUT, '>/dev/null' or die "$! opening STDOUT to /dev/null";
  open STDERR, '>/dev/null' or die "$! opening STDERR to /dev/null";

  my $bdb = $self->_open_bdb( 'tmp', 'new' );
  $bdb->initialize( -erase=>1 );

  foreach my $f ( @gff3_fh ) { #< non-verbose
    $bdb->do_load_gff($f,0);
    close $f;
  }
  foreach my $f ( @fasta_fh ) { #< non-verbose
    $bdb->load_fasta($f,0);
    close $f;
  }

  #now that we're done with bioperl, we can restore normal error reporting
  open STDERR, ">&STDERR_SAVE" or die "$! restoring STDERR";
  open STDOUT, ">&STDOUT_SAVE" or die "$! restoring STDOUT";

  #now grant web_usr select permissions
  my $bdb_dbh = $bdb->features_db();
  my $tables = $bdb_dbh->selectall_arrayref(<<EOSQL);
select tablename from pg_catalog.pg_tables where schemaname = 'public'
EOSQL
  foreach (@$tables) {
    $bdb_dbh->do(qq|grant select on "$_->[0]" to public|);
  }

  #find the new version number
  my $new_version = $self->_next_db_version();

  undef $bdb; #< close the bdb database connection so we can rename that database

  # rename the .tmp DB to the new version number
  $self->_rename_db($tmp_db,$self->dbname($new_version));

  # and now finally, make sure we don't have too many old versions
  # sitting around
  $self->_clean_up_older_databases();

  #update this object's dbname
  $self->{version} = $self->_highest_db_version();
}


#### HELPER FUNCTIONS

sub _db_exists {
  my ($self,$dbname) = @_;

  my ($exists) = $self->control_dbh->selectrow_array(<<EOSQL,undef,$dbname);
select datname from pg_catalog.pg_database where datname = ? limit 1
EOSQL

  return 1 if $exists;
  return 0;
}

sub _make_db {
  my ($self,$dbname) = @_;
  $self->control_dbh->do(<<EOSQL);
create database "$dbname"
EOSQL
}

sub _rename_db {
  my ($self, $old_name, $new_name) = @_;

  my $conns = $self->control_dbh->selectall_arrayref('select * from pg_stat_activity where datname like ?',undef,$self->dbname_root.'%');
  #print "CURRENT CONNECTIONS:\n";
  #print Dumper $conns;
  my $retry_count = my $retries = 5;
  my $success = 0;
  while($retry_count--) {
    my $result = eval {
      $self->control_dbh->do(<<EOSQL);
alter database "$old_name" rename to "$new_name"
EOSQL
    };
    unless( $@ ) {
      $success = 1;
      last;
    } else {
      warn $@;
      die unless $@ =~ /other users/;
      warn "waiting 5 minutes for $old_name to become free up...\n";
      sleep 300; #< wait 5 minutes for autovacuum to finish
    }
  }
  unless( $success ) {
    die "db rename $old_name -> $new_name failed, even after $retries tries:\n$@\n";
  }
}

#find the next version in line
sub _next_db_version {
  my ($self) = @_;
  my $d = $self->_highest_db_version();
  $d ||= 0;
  return $d+1;
}

#only keep the last couple of database versions
sub _clean_up_older_databases {
  my ($self,$num_to_keep) = @_;

  my $highest = $self->_highest_db_version
      or return;

  $num_to_keep = 3 unless defined $num_to_keep; #< by default, keep 3 db versions around

  for( my $del = $highest - $num_to_keep;
       $self->_db_exists($self->dbname_root.'.'.$del);
       $del--
     ) {

    eval { #don't care if the drop fails
      $self->control_dbh->do(qq|drop database "|.$self->dbname_root.qq|.$del"|);
    }
  }
}

sub _highest_db_version {
  my ($self) = @_;

  my $r = $self->dbname_root;
  $r =~ s!([\.\$\^])!\\$1!g;
  my ($d) = $self->control_dbh->selectrow_array(<<EOSQL,undef,$r,'^'.$r.'\.\d+$');
select regexp_replace(datname, '^' || ? || '\.','')::int as version
from pg_catalog.pg_database
where datname ~ ?
order by version desc
limit 1
EOSQL

   #warn "found highest version '$d'\n";
   return $d;
}


=head1 AUTHOR(S)

Robert Buels

=cut

__PACKAGE__->meta->make_immutable;
1;
