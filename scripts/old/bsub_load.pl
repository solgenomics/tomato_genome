#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;
use File::Spec;
use File::Temp qw/tempfile tempdir/;
use File::Basename;
use File::Copy;

use Hash::Util qw/ lock_keys lock_hash /;

use Memoize;

use Data::Dumper;

use Bio::DB::GFF;

use CXGN::DB::Connection;

use CXGN::Genomic::Clone;
use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident/;

use CXGN::Tools::List qw/ str_in /;
use CXGN::Tools::Run;
use CXGN::Tools::Script qw/lock_script unlock_script/;

use CXGN::Publish qw/parse_versioned_filepath published_as/;
use CXGN::TomatoGenome::BACPublish qw/aggregate_filename glob_pattern parse_filename/;

use CXGN::IndexedLog;

use CXGN::TomatoGenome::Config;

#########  DEFAULTS ########

my $cfg = CXGN::TomatoGenome::Config->load_locked;
our $default_log = $cfg->{'bac_loading_log'};
our $analysis_name = 'cxgnbacpipeline';
our %resource_files = (   #NOTE: these type terms have to be in SOFA
		       lycopersicum_combined_unigene_seqs => { type => 'assembly' },
		       sgn_ests_tomato_potato  => { type => 'EST' },
		       repeats_master => { type => 'repeat_family' },
		      );

# log in with the name of the user running this unless another one has been specified
$ENV{DBUSER} = '' unless defined $ENV{DBUSER};

######### /DEFAULTS ########

#note that these debugging routines return their arguments
use constant DEBUG => $ENV{BSUBLOADDEBUG} ? 1 : 0;
sub vprint(@) {
  print @_ if our $verbose or our $dry_run;
  return @_;
}
sub dprint(@) {
  print @_ if DEBUG;
  return @_;
}
sub dprinta(@) {
  dprint join(' ',@_),"\n";
  return @_;
}
sub vdo_chado(@) {
  print 'do('.join( '', map {($_||'').",\n"} @_).")\n" if our $verbose or our $dry_run;
  dbh()->do(@_) unless $dry_run;
}
sub vdo_dbgff(@) {
  print 'do('.join( '', map {($_||'').",\n"} @_).")\n" if our $verbose or our $dry_run;
  dbgff_conn()->do(@_) unless $dry_run;
}

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [ options ]  <publishing directory>

  Script to update the state of our chado database to reflect the
  state of the BAC publishing repository.  Loads sequences, their
  genbank accessions, and gff3 annotations.  To do this, it uses a log
  file (using CXGN::IndexedLog) that records when a certain 'thing'
  (file, sequence, etc) was loaded, and then figures out what it needs
  to load based on what's in the publishing repository now, and what
  the log file says has already been loaded.

  Options:

  -l  <table>
     Database table to use for logging info.
     Default: $default_log

  -x dry run.  don't actually load anything.  implies -v.

  -v be verbose

  -r reload all BAC sequences, genbank accessions, and annotations

  -R before loading BAC sequences and annotations, reinitialize the
     database by deleting all features and reloading the reference
     features (unigenes, ESTs, repeat families, etc.).  implies -r.

  -V run a vacuum analyze at the end of the load.  You can only vacuum
     tables that you own.
EOU
}

my %opt;
getopts('l:xvrRV',\%opt) or usage();
$opt{l} ||= $default_log;
my $dry_run = $opt{x};
my $verbose = $opt{v};

my $pubdir = shift @ARGV;

$pubdir && -d $pubdir or usage('publishing directory not found or not specified');

#make a lock file to keep multiple instances of this program from running
lock_script()
  or die "please do not run more than one instance of this script at the same time\n";

#open our loading log
my $dbh = CXGN::DB::Connection->new;
my $log = CXGN::IndexedLog->open( DB => $dbh, $opt{l} );

$log->reset if $opt{r} || $opt{R};

#reinitialize the DB if necessary
if( $opt{R} ) {
  $log->append("INIT DB initializing database, deleting all sequences and annotations, reloading reference sequences");
  warn "reinitializing both chado and gbrowse DBs.  This could take a really long time, please be patient...\n";
  eval {
    reinitialize_db();
  }; if( $EVAL_ERROR ) {
    $log->append("db initialization failed: $EVAL_ERROR");
    die $EVAL_ERROR;
  } else {
    $log->append("reinitialize completed");
    warn "reinitialization complete.\n";
  }
}

#load sequences and genbank accessions
#into chado
load_each( glob_pattern('all_seqs',$pubdir),
	   'load_sequence_file_chado',
	 );
#and into bio_db_gff
#load_each( glob_pattern('all_seqs',$pubdir),
#	   'load_sequence_file_bp',
#	 );

### load annotations into bio_db_gff
# #find all the *.all.gff3 files in the publishing dir,
# #and make they're all loaded
# load_each( glob_pattern('all_gff3',$pubdir),
# 	   'load_annotations',
# 	 );

#given a glob pattern and a loading subroutine, run the given
#subroutine on each of the files
sub load_each {
  my ($glob_pat,$load_fn_name) = @_;
  foreach my $file (glob $glob_pat) {
    chomp $file;

    vprint "load_each possibly loading $file with $load_fn_name...\n";

    my ($dsn) = dbh()->get_connection_parameters;
    my $spec_string = join('_',$load_fn_name,$dsn);

    #skip if this file has already been loaded
    if( $log->lookup(content => "FINISHED_LOAD_$spec_string $file") ) {
      vprint "logged as loaded, skipping.\n";
      next;
    }

    vprint "NOT logged as loaded, running $load_fn_name on it...\n";

    my $bn = basename($file);
    eval {
      no strict 'refs';
      $log->append("STARTING_LOAD_$spec_string $file") unless $dry_run;
      $load_fn_name->($file);
    }; if($EVAL_ERROR) {
      my $e = "$load_fn_name($bn) loading failed:\n$EVAL_ERROR";
      print $e;
      $log->append("FAILED_LOAD_$spec_string : $e") unless $dry_run;
      die $EVAL_ERROR if $EVAL_ERROR =~ /Got signal/;
    } else {
      $log->append("FINISHED_LOAD_$spec_string $file") unless $dry_run;
      vprint "$bn loading successful\n";
    }
  }
}

#now do a vacuum analyze if requested
if($opt{V}) {
  $log->append('running vacuum analyze') unless $dry_run;
  vprint("running vacuum analyze on the DB (this could take a LONG time)...\n");
  vdo_chado('vacuum analyze');
  vdo_dbgff('vacuum analyze');

  vprint("done.\n");
}

unlock_script(); #remove program lock file

########### SUBROUTINES ############

sub load_sequence_file_chado {
  my ($seqfile) = @_;

  my $bn = basename($seqfile);

  #check whether the sequence is in there
  #if it's not, load it
  my $p = parse_filename($seqfile)
    or die "could not parse file '$seqfile'";
  if( feature_exists($p->{seq_name}) ) {
    #one-piece seq, update its genbank accession if necessary
    update_gb_acc($p);
  }
  elsif( feature_exists( $p->{seq_name}.'-1' ) ) {
    #multi-piece seq, don't load it
  }
  else {
    ##### LOAD THE SEQ INTO CHADO
    #seq completely not in db, load it,
    #which will also load the genbank acc if present
    my $seq_gff3 = gmod_fasta2gff3($p);
    gmod_load_gff3($seq_gff3);
    #put entries into the clone_feature table for these seqs,
    #connecting them to their clone
    my $clone_id = get_clone($p->{seq_name})->clone_id;
    vdo_chado(<<EOSQL) unless $dry_run;
insert into clone_feature (clone_id,feature_id) select $clone_id,feature_id from feature f join cvterm c on(f.type_id=c.cvterm_id) where f.name like '$p->{seq_name}%' and c.name = 'BAC_clone'
EOSQL
  }
}

sub load_sequence_file_bp {
  my ($seqfile) = @_;

  my $p = parse_filename($seqfile)
    or die "could not parse file '$seqfile'";

  vprint "checking whether $p->{seq_name} exists in the Bio::DB::GFF...\n";

  ##### ALSO LOAD IT INTO bio_db_gff FOR GBROWSE
  unless( bp_feature_exists($p->{seq_name}) || bp_feature_exists($p->{seq_name}.'-1') ) {
    vprint  "sequence DOES NOT exist, loading\n";
    dbgff()->load_fasta($p->{filename})
      or die "could not load $p->{filename} into bio_db_gff database\n";
#     my $seq_gff3 = bp_fasta2gff3($p);
#     bp_load_gff3($seq_gff3);
  } else {
    vprint "sequence exists, skipping.\n";
  }
}

sub update_gb_acc {
  my ($p) = @_;
  #get the feature_id, genbank acc, and db_id and check whether there's a genbank accession
  our $genbank_db_id ||= do {
    my ($id) = dbh()->selectrow_array("select db_id from db where name = 'DB:GenBank_Accession'");
    $id or confess "no DB found with name DB:GenBank_Accession";
  };

  my ($feature_id,$genbank_acc) = dbh()->selectrow_array(<<EOSQL,undef,$genbank_db_id,$p->{seq_name});
select f.feature_id,( select dbx.accession
                      from feature_dbxref fd
                      join dbxref dbx on fd.dbxref_id=dbx.dbxref_id
                      where dbx.db_id = ?
                        and fd.feature_id = f.feature_id
                    )
from feature f
where f.name = ?
EOSQL

  unless($genbank_acc) {
    #insert a genbank accession dbxref
    my $pdl = parse_seq_file_defline($p->{filename});

    return unless $pdl->{gb_accession} && $pdl->{gb_version};

    #make a dbxref
    vdo_chado(<<EOSQL,undef,$genbank_db_id,$pdl->{gb_accession},$pdl->{gb_version});
insert into dbxref
       (db_id,accession,version)
values (?,?,?)
EOSQL

    #link it to our feature
    vdo_chado(<<EOSQL,undef,$feature_id,dbh()->last_insert_id('dbxref'));
insert into feature_dbxref
       ( feature_id, dbxref_id )
values ( ?,          ?         )
EOSQL
  }
}

#   #if it is, just update its genbank accession
#   my $seq_in = Bio::SeqIO->new(-file   => $seqfile,
# 			       -format => 'fasta',
# 			      );

#   while( my $seq = $seq_in->next_seq ) {
#     my $clone = get_clone($seq->primary_id);
#     my ($acc) = split /\s+/, $seq->desc;
#     $acc =~ /^[A-Z_]{2,3}\d+\.\d+$/ or die "'$acc' does not look like a valid genbank accession\n";

#     our $tomato_organism_id ||= do{
#       my ($id) = dbh()->selectrow_array("select organism_id from organism where common_name='tomato'");
#       $id
#     };
#     our $bac_type_id ||= do {
#       my ($id) = dbh()->selectrow_array("select cvterm_id from cvterm where name='BAC_clone'");
#       $id
#     };

#     #insert the sequence if necessary
#     unless( dbh()->selectrow_arrayref("select feature_id from public.feature where name=? and type_id=?",
# 				      undef,
# 				      $seq->primary_id,
# 				      $bac_type_id
# 				     )
# 	  ) {
#       #load the sequence into the feature table
#       dbh()->do(<<EOSQL,undef, $tomato_organism_id, ($seq->primary_id) x 2, $seq->seq, $seq->length, $bac_type_id);
# insert into public.feature
#        (organism_id, name, uniquename, residues, seqlen, type_id)
# values (?,           ?,    ?,          ?,        ?,      ?      )
# EOSQL
#     }

#   }

# }

memoize('get_clone');
sub get_clone {
  CXGN::Genomic::Clone->retrieve_from_clone_name(@_)
      or die "could not retrieve clone from name @_";
}


sub load_annotations {
  my ($gff3file) = @_;

  my $bn = basename($gff3file);

  vprint "loading $bn\n";

  #parse all the BAC and version information from the filename
  my $pub = CXGN::TomatoGenome::BACPublish::publisher();
  my $vparsed = $pub->parse_versioned_filepath($gff3file)
    or die "Could not parse versioned filename '$gff3file'";
  my $cparsed = parse_clone_ident($vparsed->{name},'versioned_bac_seq')
    or die "Could not parse BAC name '$vparsed->{name}'";
  my $clone = CXGN::Genomic::Clone->retrieve_from_parsed_name($cparsed)
    or die "Could not find BAC in database with name '$vparsed->{name}'";

  #is this bac sequence already in the chado db?
  die "Aborting load, BAC sequence $vparsed->{name} not in chado DB (does a new sequence version exist?)"
    unless ( feature_exists($vparsed->{name}) || feature_exists("$vparsed->{name}-1") );

  #is this bac sequence already in the bioperl db?
  unless ( bp_feature_exists($vparsed->{name}) || bp_feature_exists("$vparsed->{name}-1") ) {
    if( feature_exists($vparsed->{name}) ) {
      vprint "transferring $vparsed->{name} sequence from chado into bio_db_gff...";
      my ($seq) = dbh()->selectrow_array('select residues from feature where name = ?',undef,$vparsed->{name})
	or die "error.  No sequence in chado for $vparsed->{name}.\n";
      dbgff()->load_sequence_string($vparsed->{name},$seq);
      vprint "done.\n";
    }
    else {
      die "Aborting load, BAC sequence $vparsed->{name} not in bio_db_gff DB (does a new sequence version exist?)"
    }
  }


  #delete any annotations it has
  delete_features($vparsed->{name});

  #now finally load these annotations
  #  gmod_load_gff3($gff3file, '--analysis' => $analysis_name);
  bp_load_gff3($gff3file); #< only load annotations into the gbrowse database
}

# given a seq file, return a hashref with info parsed out of it,
# with keys gb_accession, gb_version, htgs_phase, and sequenced_by,
# not all of which are necessarily present
sub parse_seq_file_defline {
  my ($seqfile) = @_;
  my $seq = Bio::SeqIO->new(-file => $seqfile, -format => 'fasta')->next_seq;

  lock_keys( my %info, qw/ htgs_phase gb_accession gb_version sequenced_by upload_account_name /);

  if($seq->desc =~ /htgs_phase:(\d)/) {
    $info{htgs_phase} = $1;
    die "invalid phase $info{htgs_phase}" unless str_in($info{htgs_phase},1..3);
  }

  if($seq->desc =~ /^[A-Z_]{2,3}\d+\.(\d+)/) {
    $info{gb_accession} = $MATCH;
    $info{gb_version}   = $1;
  }

  if( $seq->desc =~ /sequenced_by:(\S+)/ ) {
    $info{sequenced_by} = $1;
  }

  if( $seq->desc =~ /upload_account_name:(\S+)/ ) {
    $info{upload_account_name} = $1;
  }

  #lock_hash(%info);

  return \%info;
}

#args: a parsed filename, optional hash-style list of options
#ret:  a gff3 filename with those seq
#converts the given fasta into gff3, in a temp file
sub gmod_fasta2gff3 {
  my ($p,%args) = @_;

  my $bn = $p->{basename};

  vprint "converting $bn to gmod gff3...\n";

  #make a temp dir to hold all this crap
  my $tempdir = tempdir( File::Spec->catdir(File::Spec->tmpdir,
					    'bsub-load-bac-annotations-gmod_fasta2gff3-XXXXXXXX',
					   ),
			 CLEANUP => !DEBUG);

  #make a subdir to hold our temp seq file, symlink our seq file there
  my $seq_dir = File::Spec->catdir($tempdir,'tempseq');
  mkdir $seq_dir;
  my $temp_seq_file = File::Spec->catfile($seq_dir,'temp.fasta');
  system("ln -s $p->{filename} $temp_seq_file");
  -r $temp_seq_file or die "no temp seq file '$temp_seq_file'";

  #figure out some facts about the seq file
  my $pdl = parse_seq_file_defline($temp_seq_file);
  #my ($acc,$ver,$phase) = parse_seq_file_defline($temp_seq_file);

  #now make a temp file to hold the gff3 output
  my $gff3_file = File::Spec->catfile($tempdir,$bn);

  my $attr = join( ';',
		   ( $pdl->{gb_accession}        ? "Dbxref=GenBank_Accession:$pdl->{gb_accession}"      : () ),
		   ( defined($p->{finished})     ? ("finished_seq=$p->{finished}")                      : () ),
		   ( $pdl->{htgs_phase}          ? ("htgs_phase=$pdl->{htgs_phase}")                    : () ),
		   ( $pdl->{sequenced_by}        ? ("sequenced_by=$pdl->{sequenced_by}")                : () ),
		   ( $pdl->{upload_account_name} ? ("upload_account_name=$pdl->{upload_account_name}")  : () ),
		 );

  CXGN::Tools::Run->run('gmod_fasta2gff3.pl',
			'--fasta_dir'   => $seq_dir,
			'--gfffilename' => $gff3_file,
			'--type'        => $args{type}   || 'BAC_clone',
			'--source'      => $args{source} || 'SGN',
			$attr ? ('--attributes'  => $attr) : (),
		       );

  -s $gff3_file or die "$bn fasta->gff3 conversion produced a gff3 file of zero size";
  vprint "done.\n";
  return $gff3_file;
}

# #args: a parsed filename, optional hash-style list of options
# #ret:  a gff3 filename with those seq
# #converts the given fasta into gff3, in a temp file
# sub bp_fasta2gff3 {
#   my ($p,%args) = @_;

#   my $bn = $p->{basename};

#   vprint "converting $bn to bp gff3...\n";

#   my $seqlen = Bio::SeqIO->new(-file => $p->{filename}, -format => 'fasta')->next_seq->length;

#   #now make a temp file to hold the gff3 output
#   my ($gff3fh,$gff3_file) = tempfile(File::Spec->tmpdir.'/bp_fasta2gff3-XXXXXXXX', UNLINK => !DEBUG)
#     or die "cannot make tempfile for bp_fasta2gff3";

#   open my $seqfh,$p->{filename} or die "$! opening $p->{filename}\n";

#   print $gff3fh "##gff-version 3\n##sequence-region $p->{seq_name} 1 $seqlen\n";
#   print $gff3fh $_ while <$seqfh>;
#   close $gff3fh;

#   -s $gff3_file or die "$bn fasta->bp gff3 conversion produced a gff3 file of zero size";
#   vprint "done.\n";
#   return $gff3_file;
# }

#args: a gff3 file
#loads a gff3 file into our chado db
sub gmod_load_gff3 {
    my ($file,%args) = @_;
    my $bn = basename($file);
    vprint "loading gff3 file $bn into chado database...\n";
    $file or die 'must pass a file to load';
    unless( our $dry_run ) {
        CXGN::Tools::Run->run('gmod_bulk_load_gff3.pl',
                              '--gfffile'  => $file,
                              '--organism' => 'tomato',
                              '--dbname' => dbh()->dbname,
                              '--dbuser' => dbh()->dbuser, #< do not pass user and password to gff3 loader
                              '--dbpass' => dbh()->dbpass,
                              '--dbhost' => dbh()->dbhost,
                              '--dbport' => dbh()->dbport,
                              '--skip_vacuum', #< we'll do a vacuum analyze when we're done with all the loads
			      '--recreate_cache',
			      ( DEBUG ? '--save_tmpfiles' : () ),
                              %args,
                              {
                               working_dir => File::Spec->tmpdir }
                             );
    } else {
        vprint "dry run, skipping load of $bn\n";
    }
    vprint "done.\n";
}

#args: a gff3 file
#loads a gff3 file into our bio::db::gff db
sub bp_load_gff3 {
  my ($file) = @_;
  my $bn = basename($file);
  vprint "loading gff3 file $bn into Bio::DB::GFF database...\n";
  $file or die 'must pass a file to load';

  #filter the GFF3 to qualify the class of the Targets if they are bac sequences
  my ($gff3fh,$gff3_file) = tempfile(File::Spec->tmpdir.'/bp_load_gff3-XXXXXXXX', UNLINK => !DEBUG)
    or die "cannot make tempfile for bp_load_gff3";
  open my $in_fh,$file or die "$! opening $file";
  while (<$in_fh>) {
    if (/Target=([^\s;]+)/ && parse_clone_ident($1,'versioned_bac_seq') ) {
      s/Target=/Target=Sequence:/;
    }
    print $gff3fh $_;
  }
  close $gff3fh;
  close $in_fh;

  unless( our $dry_run ) {
    my ($dsn, $user, $pass) = dbgff_conn_params();

    CXGN::Tools::Run->run('bp_load_gff.pl',
			  -d          => $dsn,
			  '--user'    => $user,
			  '--pass'    => $pass,
			  '--adaptor' => 'dbi::pg',
			  $gff3_file
			 );
  } else {
    vprint "dry run, skipping bp load of $bn\n";
  }
  vprint "done.\n";
}

#args: a feature name
#delete all the subfeatures of the given feature from the chado db
#but don't delete the actual parent feature
sub delete_features {
  my ($ref_seq_name) = @_;

  my $clone = CXGN::Genomic::Clone->retrieve_from_clone_name($ref_seq_name)
    or die "could not retrieve clone for seq name '$ref_seq_name'";

  #delete subfeatures
  vprint "deleting subfeatures for '$ref_seq_name'...";
  vdo_chado(<<EOQ,undef,$ref_seq_name,'^'.$ref_seq_name.'-\d+'."\$");
delete from feature
where feature_id in(
  select feature_id
  from featureloc
  where srcfeature_id in(
    select feature_id from feature where name=? or name ~ ?
  )
);
EOQ

  delete_chado_load_cache(); #need to do this whenever we change the features;

  my $dbgff = dbgff();
  unless($dry_run) {
    #save the reference sequence in memory
    my $segment = $dbgff->segment($ref_seq_name);
    my $seq;
    $seq = $segment->seq if $segment;

    $dbgff->delete(-ref => $ref_seq_name);
    $dbgff->delete(-ref => $_) foreach @{dbh()->selectcol_arrayref(<<EOQ,undef,$ref_seq_name.'-\d+'."\$")};
select name from feature where name ~ ?
EOQ
    $dbgff->load_sequence_string($seq->display_id,$seq->seq) if $seq; #< put the reference sequence back
  }

  vprint "done.\n";
}

#cached database handle, created the first time it's needed
sub dbh {
  return our $_dbconn ||= CXGN::DB::Connection->new({ config => $cfg, dbargs => {AutoCommit=>1} });
}

sub dbgff_conn_params {
  my ($dsn, $user, $pass) = dbh()->get_connection_parameters;
  $dsn =~ s/dbname=\w+/dbname=bio_db_gff/;
  return ($dsn,$user,$pass);
}

sub dbgff_conn {
  return our $_dbgff_conn ||= do {
    my ($dsn, $user, $pass) = dbgff_conn_params();
    DBI->connect($dsn,$user,$pass,{AutoCommit => 1})
      or die "could not connect to Bio::DB::GFF db using dsn '$dsn'"
  }
}

sub dbgff {
  return our $_dbgff ||= do {
    my ($dsn, $user, $pass) = dbgff_conn_params();
#    warn "connecting with $dsn,$user,$pass\n";
    Bio::DB::GFF->new( -adaptor => 'dbi::pg',
		       -dsn => $dsn,
		       -user => $user,
		       -pass => $pass,
		     )
  }
}


sub reset_sequences {
  my %sequences = (
		   feature              => "feature_feature_id_seq",
		   featureloc           => "featureloc_featureloc_id_seq",
		   feature_relationship => "feature_relationship_feature_relationship_id_seq",
		   featureprop          => "featureprop_featureprop_id_seq",
		   feature_cvterm       => "feature_cvterm_feature_cvterm_id_seq",
		   dbxref               => "dbxref_dbxref_id_seq",
		   synonym              => "synonym_synonym_id_seq",
		   feature_synonym      => "feature_synonym_feature_synonym_id_seq",
		   feature_dbxref       => "feature_dbxref_feature_dbxref_id_seq",
		   analysisfeature      => "analysisfeature_analysisfeature_id_seq"
		  );

  while (my ($table,$seqname) = each %sequences) {
    vprint "resetting sequence for $table table...\n";
    vdo_chado("select setval('$seqname',(select max(${table}_id)+1 from $table),false)");
  }
}

#args: none
#delete all features from the database and load the source features that are
#used in our analyses
sub reinitialize_db {

  vprint("reinitializing database, deleting everything from '$analysis_name' analysis\n");

  vdo_chado(<<EOQ,undef,$analysis_name);
delete from feature
where feature_id in( select feature_id
                     from analysisfeature
                     join analysis using(analysis_id)
                     where name=?
                   )
EOQ

  #delete all orphan synonyms
  vdo_chado(<<EOQ,undef,$analysis_name);
delete from synonym
where synonym_id in( select synonym_id
                     from synonym
                     left join feature_synonym using(synonym_id)
                     where feature_synonym_id is null
                   )
EOQ

  vdo_chado("delete from analysis where name=?",undef,$analysis_name);

  reset_sequences();
  delete_chado_load_cache();

  vdo_chado(<<EOQ,undef,$analysis_name,$analysis_name);
insert into analysis
(name,description,program,programversion)
values
(?,'the SGN BAC analysis pipeline',?,'any')
EOQ

  #delete everything from the bio_db_gff db
  unless( our $dry_run ) {
    vprint "deleting all features from bio_db_gff db...";
    dbgff()->delete(-force => 1);
    vdo_dbgff('delete from fdna'); #the above delete does not clean up fdna
    vprint "done.\n";
  } else {
    vprint "dry run, Bio::DB::GFF recreate\n";
  }

  #now load all of our resource files into chado
  while( my ($filetag, $args) = each %resource_files ) {
    vprint "fetching resource file $filetag...\n";
    unless($dry_run) {
      my $local_resource_file = CXGN::TomatoGenome::BACPublish::resource_file($filetag);
      my ($bn,$dir) = fileparse($local_resource_file);
      my $p = { basename => $bn, dirname => $dir, filename => $local_resource_file };
      my $gff3_resource = gmod_fasta2gff3($p, %$args);
      gmod_load_gff3($gff3_resource, '--analysis' => $analysis_name);
    }
    vprint "done.\n";
  }
}


=head2 feature_exists

  Usage: feature_exists('some_feature_name');
  Desc : look to see whether a feature with the given name exists in our chado db
  Ret  : the feature's name if it is there, undef otherwise
  Args : a feature name, (optional) dbh to use
  Side Effects:  looks up things in the chado db
  Example:

=cut

sub feature_exists {
  my ($feature_name) = @_;

  my $dbh = dbh();

  my ($cnt) = $dbh->selectrow_array('select count(*) from public.feature where name=?',undef,$feature_name);
  return $cnt > 0 ? 1 : undef;
}

sub bp_feature_exists {
  my ($featurename) = @_;
  my $seg = dbgff->segment($featurename);
  return $seg ? 1 : 0;
}

#chado keeps a materialized view of the feature table
sub delete_chado_load_cache {
  dbh()->dbh_param(PrintError => 0);
  eval { vdo_chado("drop table tmp_gff_load_cache"); }; #don't care if that fails
  dbh()->dbh_param(PrintError => 1);
}
