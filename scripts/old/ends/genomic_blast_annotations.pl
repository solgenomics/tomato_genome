#!/usr/bin/perl
use strict;
use warnings;
use English;

#standard perl and cpan
use Getopt::Std;
use File::Temp qw/ tempdir /;
use File::Basename;
use File::Spec;
use FindBin;

#vanilla BioPerl
use Bio::Tools::Run::StandAloneBlast;
use Bio::SearchIO;
use Bio::SearchIO::FastHitEventBuilder;

#DB connection
use CXGN::DB::Connection;

#CXGN configuration system
use CXGN::TomatoGenome::Config;

#Genomic framework
use CXGN::Genomic::Library;
use CXGN::Genomic::GSS;
use CXGN::Genomic::Search::GSS;

#custom CXGN BioPerl
use Bio::SeqIO::xform;
use Bio::Tools::Run::CXGNRemoteBlast;

#other
use CXGN::Annotation::BlastAnnotator;

CXGN::DB::Connection->verbose(0); #turn off annoying DB::Connection messages
#CXGN::CDBI::Class::DBI->db_Main->trace(2,'/tmp/dbitrace');
CXGN::CDBI::Class::DBI->db_Main->dbh_param( AutoCommit => 1 ); #using localized transactions


########## CONFIGURATION VARS ###############
my $cxgn_conf = CXGN::TomatoGenome::Config->new;

my $default_databases_root = $cxgn_conf->get_conf('blast_db_path') || '/data/shared/blast/databases/current';
my $default_temp_directory = '/data/local';

my $blast_evalue = 1e-10;

my $parallel_blast_host = 'amatxu.sgn.cornell.edu';

#maximum number of hits to store in the database
#for a single BLAST query
my $max_stored_hits = 10;

#minimum size of a BLAST DB (in sequences) for a remote BLAST to be executed
my $min_remote_db_size    = 20_000;
my $min_remote_query_size = 20_000;

#############################################


sub usage {
  die <<EIEIO;

Usage:

  $FindBin::Script [ options ]

  -t   Directory in which to create (and clean up) temporary working directory.
       Defaults to $default_temp_directory

  -d   <blast_database_path>
       Base path to BLAST databases. Default $default_databases_root.

  -r   Allow remote parallel blasts on the CXGN computation cluster.

  -L   <shortname>,<shortname>
       Only run annotations against libraries with the given shortnames.
       By default, runs any necessary annotations and contam screens for all
       Genomic libraries.

  -e   Kludgily use one or more externally-generated BLAST reports.
       See below for an explanation.

  -c   Annotate with only contamination BLAST databases, skipping
       the (large) annotation databases (like nr)

  -b <blast db name>,<blast db name>
       Only annotate against the blast DBs whose file_base's appear
       in this list.

  USING EXTERNALLY-GENERATED BLAST REPORTS
    If the -e option is specified, this script will create its temporary
    working directory, tell you where that is, and tell you
    exactly where you can put (or symlink) your BLAST reports for it to use,
    if you have any.

    NOTE:  The genomic.blast_query.source_id field is taken directly from the
           name of each query sequence in the blast report.  Make sure that
           the query names in your blast report are numerical and correspond
           to the GSS IDs of the query sequences, or else this script WILL NOT
           DO THE RIGHT THING.

EIEIO
}
#store command-line options in the %opt hash
our %opt;
getopts('ed:t:rcb:L:',\%opt) or usage();

sub progprint(@); #unbuffered progress message.  defined below.



if($opt{r} && !$opt{t}) {
  #they requested to run blasts on the cluster, but did not specify a temp directory
  print "WARNING: '-r' specified without '-t'.  Are you sure $default_temp_directory is readable from the CXGN cluster?\n";
}
if($opt{r} && !$opt{d}) {
  #they requested to run blasts on the cluster, but did not specify a directory for BLAST databases
  print "WARNING: '-r' specified without '-d'.  Are you sure $default_databases_root is readable from the CXGN cluster?\n";
}

my $blast_databases_root = $opt{d} ? File::Spec->rel2abs($opt{d}) : $default_databases_root;

### connect to the development database ###
my $dbconn = CXGN::DB::Connection->new('genomic');

### make a temporary dir in the working directory ###
$opt{t} ||= $default_temp_directory; #default temp directory
-w $opt{t}
  or die "Cannot write to temp directory $opt{t}\n";
$opt{t} = File::Spec->rel2abs($opt{t});

my $tempdir = tempdir( File::Spec->catdir( $opt{t}, 'genomic-blast-annotations-XXXXXXX'),
		       CLEANUP => 0 )
  or die "Could not create temporary directory in $opt{t} ($!)\n";


### instantiate a BlastAnnotator ###
my $annotator = CXGN::Annotation::BlastAnnotator->new;
$annotator->source_type_shortname('gss');


### assemble the list of genomic libraries to update the blast annotations on ###
my @libs = do {
  if( $opt{L} ) {
    map { CXGN::Genomic::Library->search( shortname => $_ ) } (split ',', $opt{L});
  }
  else {
    CXGN::Genomic::Library->retrieve_all
  }
};

@libs or die "No libraries found.\n";

if( $opt{e} ) {
  external_reports_prompt( $tempdir, $blast_databases_root, @libs );
}

#for each library in Genomic
foreach my $lib (@libs) {

  progprint "Annotating library ".$lib->shortname."\n";

  #look up the blast databases for this library
  my @blastdbs = get_blast_dbs_for_lib($lib,$blast_databases_root);

  #   for each blast DB they need to be screened/annotated against
  foreach my $bdb_record (@blastdbs) {
    my ($bdb,$is_contaminant) = @$bdb_record;

    # get a Bio::SearchIO parser handle, connected toeither an externally-generated blast report
    # or the results of a new BLAST run
    my $searchio_parser = do {
      if($opt{e} && -f (my $ext_report_name = blast_report_name($tempdir,$lib,$bdb)) ) {
	progprint "   using external BLAST report $ext_report_name\n";
	parse_ext_report($ext_report_name);
      } else {
	progprint "   with BLAST database ".$bdb->subdir.'/'.$bdb->file_basename."\n";
	progprint blast_report_name($tempdir,$lib,$bdb)." not provided, doing a new blast...\n" if $opt{e};
	run_search($tempdir,$lib,$bdb_record);
      }
    };

    if( UNIVERSAL::isa($searchio_parser,'Bio::SearchIO') ) {
      eval { #record the hits for this report inside a transaction
	CXGN::CDBI::Class::DBI->db_Main->begin_work;
	$annotator->record_hits( $bdb, $searchio_parser,
				 $is_contaminant ? (\&flag_as_contaminated, \&flag_as_contam_screened) : undef,
			       );
	#print "Type 'yes' to commit this transaction... ";
	#die "Commit aborted\n" unless <STDIN> =~ /^yes/i;
	CXGN::CDBI::Class::DBI->db_Main->commit;
	progprint "Committed blast report load transaction.\n";
      }; if( $EVAL_ERROR ) {
	my $error = $EVAL_ERROR;
	eval { #the rollback could fail also!
	  warn "ERROR: Failed to record hits for BlastDB ".$bdb->file_basename.": $error\n";
	  CXGN::CDBI::Class::DBI->db_Main->rollback;
	}; if( $EVAL_ERROR ) {
	  die "FAILED to rollback transaction while loading blast report!";
	} else {
	  warn "Successfully rolled back failed load.\n";
	}
      }
    }
  } #each blast DB
}#each library

#if a hit is seen on a contaminant database, flag that GSS as contaminated
our %already_flagged_as_contaminated;
sub flag_as_contaminated {
  my ($result,$hit) = @_;

  my $gssid = $result->query_name;
  unless( $already_flagged_as_contaminated{$gssid} ) {
    #avoid the overhead of flagging something multiple times
#    warn "Flagging GSS $gssid as contaminated\n";

    my $gss = CXGN::Genomic::GSS->retrieve($gssid)
      or die "No GSS object found with id $gssid";
    $gss->set_flags('hostcontam');
    $gss->update;
  }
}

#if a GSS is blasted against a contaminant database, set its status as having been screened
#for contaminants
sub flag_as_contam_screened {
  my $result = shift;
  my $gssid = $result->query_name;
  #warn "Flagging GSS $gssid as contam screened\n";

  my $gss = CXGN::Genomic::GSS->retrieve($gssid)
    or die "No GSS object found with id $gssid";
  $gss->unset_status('contam_unk');
  $gss->update;
}


{ #cached sourcetype
  my $gss_stype;
  sub gss_sourcetype {
    $gss_stype ||= shift->get_sourcetype_id_by_shortname('gss');
  }
}

#unbuffered progress message print
sub progprint(@) {
  local $| = 1;
  print @_;
}

=head2 blastprog_to_file_ext

  Usage: my $ext = blastprog_to_file_ext( 'blastn' );
         #will return '.nsq'
  Desc : get the proper file extension for one of the component
         files in a b
  Ret  :
  Args :
  Side Effects:
  Example:

=cut

sub blastprog_to_file_ext {
  my $blastprog = shift;
  if( $blastprog =~ /blastn|tblastn|tblastx/ ) {
    return '.nsq';
  } elsif( $blastprog =~ /blastp|blastx/ ) {
    return '.psq';
  }
  die "Invalid BLAST program '$blastprog'";
}

=head2 sequences_in_blast_database

  Args: the full basename of a BLAST database
  Ret : the number of sequences in the blast database with the given basename

=cut

sub sequences_in_blast_database {
  my ($basename) = @_;

  my $stats = `fastacmd -d $basename -I`;
  #  print "Got stats $stats\n";

  my ($numseqs) = $stats =~ /((?:\d+,)*\d+) sequences;/;
  $numseqs ||= 0;
  $numseqs =~ s/,//g;
  return $numseqs;
}


=head2 run_search

  Desc: run a new SearchIO (BLAST) search of the seqs in a given library
        against a given Blast database
  Args: CXGN::Genomic::Library,  CXGN::Genomic::BlastDB object to search
  Ret : BioPerl SearchIO object to iterate over the resulting SearchIO report

=cut

sub run_search {
  my ($tempdir,$lib,$bdb_record) = @_;

  my ($bdb,$is_contaminant,$dbseqs,$db_full_basename) = @$bdb_record;

  my $libid = $lib->library_id;

  ### modification time of this database ###
  my $blastprog = $bdb->blast_program;
  my $type_ext = blastprog_to_file_ext($blastprog);
  my $db_mtime = (stat($db_full_basename.$type_ext))[9];

  die "BLAST database file not found: $db_full_basename$type_ext\n"
    unless $db_mtime && (-f $db_full_basename.$type_ext);

  # find sequences that need to be done/redone against this database
  my $gss_search = CXGN::Genomic::Search::GSS->new;

  progprint "    Searching for GSSs to annotate against ".$bdb->file_basename."\n";
  my $gss_query = $gss_search->new_query;
  #  $gss_query->gss_id(" < 300"); #test, for speed
  # find sequences only in this library
  $gss_query->library_id("= $libid");
  $gss_query->flags_not_set('complexity | error | anomaly');
  $gss_query->status_not_set('legacy | discarded | deprecated | censored | vec_unk');

  # find sequences that either have no query against this DB or
  # whose query is older than the DB timestamp
  $gss_query->needs_auto_annot($bdb->blast_db_id,$db_mtime);

  #      BLAST them, BLAST them I say!
  my $seqs_file = File::Spec->catfile($tempdir,
				      'blast_query_'.$lib->shortname.'_vs_'.$bdb->file_basename.'.seq');
  my $queryseqs = fast_dump_gss_seqs($gss_search,$gss_query,$seqs_file);

  unless($queryseqs > 0) { #stop now if there are no sequences to blast
    progprint "     no sequences need annotation, skipping BLAST against ".$bdb->file_basename."...\n";
    return;
  }

  progprint "      running BLAST against ".$bdb->file_basename." ($dbseqs sequences)...\n";

  my $reportname = blast_report_name($tempdir,$lib,$bdb);

  if ( !$opt{r} && $dbseqs >= $min_remote_db_size ) {
      print "WARNING: ".$bdb->file_basename." is big ($dbseqs seqs), and you haven't allowed remote alignment searches (blasts) with the '-r' option.  You might want to consider doing that.\n";
  } elsif (!$opt{r} && $queryseqs >= $min_remote_query_size) {
      print "WARNING: This set of input sequences is big ($queryseqs seqs), and you haven't allowed remote alignment searches (blasts) with the '-r' option.  You might want to consider doing that.\n";
  }

  my $blastfactory;
  if ($dbseqs >= $min_remote_db_size) {
      print "      running remote blast ($parallel_blast_host)\n";
      $blastfactory = Bio::Tools::Run::CXGNRemoteBlast->new( database => $db_full_basename,
							     outfile  => $reportname,
							     program  => $blastprog,
							   );
      $blastfactory->remotehost($parallel_blast_host);
  } elsif ($opt{r} && $queryseqs >= $min_remote_query_size) {
      print "      running remote blast ($parallel_blast_host)\n";
    $blastfactory = Bio::Tools::Run::CXGNRemoteBlast->new( database => $db_full_basename,
							   outfile  => $reportname,
							   program  => $blastprog,
							 );
    $blastfactory->remotehost($parallel_blast_host);
  } else {
    print "     running local blast\n";
    $blastfactory = Bio::Tools::Run::StandAloneBlast->new( database => $db_full_basename,
							   outfile  => $reportname,
							   program  => $blastprog,
							 );
  }
#  $blastfactory->verbose(1);

  $blastfactory->e($blast_evalue);
  $blastfactory->v(100); #limit to best 100 hits
  $blastfactory->b(100); #limit to best 100 hits
  return $blastfactory->blastall($seqs_file);
}

=head2 seqio_dump_gss_seqs

  Usage:
  Desc :
  Ret  : number of sequences written to seq file
  Args : gss search object, gss query object, filename to dump to
  Side Effects: opens the file, dumps the sequences into it
  Example:

=cut

sub seqio_dump_gss_seqs {
  my ($gss_search, $gss_query, $seqtemp) = @_;

  $gss_search->page_size(60_000);
  my $gss_seqio = $gss_search->seqIO_search($gss_query);

  #Bio::SeqIO::xform is to SeqIOs what the perl 'map' function is to arrays/lists
  my $trimmed_seqio =
    Bio::SeqIO::xform->new(-enclose  => $gss_seqio,
			   -map_next => sub {
			     my $orig = shift;
			     my $s = $orig->vector_trimmed_trunc;
			     $s->display_id($s->primary_id);
			     $s->desc('');
			     $s
			   });

  #count the number of sequences in the database and in the query set,
  #so we can make decisions about how best to do the BLAST
  my $queryseqs = $gss_seqio->result->total_results;

  my $out = Bio::SeqIO->new(-file   => ">$seqtemp",
			    -format => 'fasta',
			    -flush  => 0,
			   );
  progprint "      writing $queryseqs query sequences to $seqtemp ...\n";
  while ( my $seq = $trimmed_seqio->next_seq ) {
    $out->write_seq($seq);
  }
  progprint "done\n";

  return $queryseqs;
}

=head2 fast_dump_gss_seqs

Same as seqio_dump_gss_seqs, but written in a much speedier way.

=cut

sub fast_dump_gss_seqs {
  my ($gss_search, $gss_query, $seqtemp) = @_;
  $gss_search->page_size(60_000);
  my $gss_results = $gss_search->do_search($gss_query);
  $gss_results->autopage($gss_query,$gss_search);

  my $queryseqs = $gss_results->total_results;
  progprint "      writing $queryseqs query sequences to $seqtemp ...\n";

  open( my $fh, ">$seqtemp" )
    or die "Could not open '$seqtemp' for writing\n";
  while( my $gss = $gss_results->next_result ) {
    print $fh '>',$gss->gss_id,"\n",$gss->trimmed_seq,"\n";
  }

  return $queryseqs;
}


=head2 parse_ext_report

  Desc:	
  Args:	full path to the blast report file
  Ret :	Bio::SearchIO parser object to parse the given blast report

=cut

sub parse_ext_report {
  my ($reportname) = @_;

  Bio::SearchIO->new( -format => 'blast',
		      -file   => $reportname );
}

=head2 blast_report_name

  Desc:	
  Args:	containing directory, library object, blastDB object
  Ret :	proper filename to use for a blast report

=cut

sub blast_report_name {
  my ($dirname, $lib, $bdb) = @_;
  return File::Spec->catfile( $dirname, 'blast_report_'.$lib->shortname.'_vs_'.$bdb->file_basename.'.blast');
}


=head2 external_reports_prompt

  Desc:	print a prompt for the user to copy in any external BLAST
        reports they may have
  Args:	temporary directory where the external reports should be put,
        a restore()'d Library object to use for iterating over the
        libraries
  Ret :	nothing useful

=cut

sub external_reports_prompt {
  my ($tempdir,$blast_databases_root,@libs) = @_;

  my $reportnames;
  foreach my $lib (@libs) {
    #look up the blast databases for this library
    my @blastdbs = map { $_->[0] } get_blast_dbs_for_lib($lib,$blast_databases_root);

    foreach my $bdb (@blastdbs) {
      $reportnames .= blast_report_name($tempdir,$lib,$bdb)."\n";
    }
  }

  progprint <<EOT;

External reports option selected.  If you have any BLAST reports for these libraries vs. these databases, you may copy them to any of the following:
$reportnames

After you are finished copying, press enter to continue.
EOT
  my $throwaway = <STDIN>
}


=head2 get_blast_dbs_for_lib

  Usage: my @dbs = get_blast_dbs_for_lib($lib,'/data/shared/blast/databases/current');
  Desc : get a list of the Blast DBs to annotate a given library against
  Ret  : list as:
         ( [ BlastDB object,
             1 if contaminant,
             number of sequences in blast database,
             full basename of blast database files,
           ],
           ...
         )
  Args : library name, name of blast databases root dir
  Side Effects: SELECTs from the database in various ways
  Example:

=cut

sub get_blast_dbs_for_lib {
  my ($lib,$blast_databases_root) = @_;

  my @blastdbs;
  if ( $opt{b} ) {
    my @names = split ',',$opt{b}
      or usage;
    @blastdbs = map { CXGN::BlastDB->search( file_base => $_ ) } @names;
    @blastdbs = map { [ $_, $_->is_contaminant_for($lib) ] } @blastdbs;
  } else {
    @blastdbs = ((map {[$_,1]} $lib->contamination_blast_dbs()),
		 #only add annotation dbs if we didn't get a -c
		 $opt{c} ? () : (map {[$_,0]} $lib->annotation_blast_dbs()),
		);
  }

  #calculate the full pathname and size in seqs of each database
  #and add them to the array
  @blastdbs = map {
    my $db_full_basename = File::Spec->catfile( $blast_databases_root,
						$_->[0]->subdir,
						$_->[0]->file_basename
					      );
    my $dbseqs = sequences_in_blast_database($db_full_basename);
    [@$_, $dbseqs, $db_full_basename]
  } @blastdbs;

  #sort our blast databases by their size (in seqs), so we blast against the smallest ones first
  @blastdbs = sort {$a->[2] <=> $b->[2]} @blastdbs;

  return @blastdbs;
}
