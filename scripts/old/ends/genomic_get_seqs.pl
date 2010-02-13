#!/usr/bin/perl
use strict;
use warnings;

use FindBin;

use English;
use DBI;
use Getopt::Std;

use CXGN::DB::Connection;
use CXGN::Genomic::GSS;

######## CONFIGURATION VARS ########

$|=1; #don't buffer progress text

#no other configuration vars right now

####################################

### print usage info and exit ###
sub usage {
   die <<EIEIO;

  Program to query the genomic database and produce two FASTA files, one of
  of sequences and one of quality values, with the given basename.

  Examples:
     $FindBin::Script -g foo

       will produce a FASTA-format file foo.seq in the current directory
       containing all genomic sequences in the database, with GSS IDs as
       identifiers.

     $FindBin::Script -g -q foo

       will produce the same .seq file, but will also produce a corresponding
       foo.qual file

  NOTE: Due to a bug in the current version of BioPerl, any .qual files produced
        will be missing the comment portions of the identifier lines.


  Usage: query-genomic-seqs.pl [ options ] <output basename>

  -T        Do not trim sequences and quals for vector and quality (output raw seqs).
  -g        Output only GSS ids as identifiers, with no comments.
            This option is useful when querying sequences for blast_vector_lexer.pl
	    or contaminant screens.
  -l <len>  Only output sequences with at least <len> bases.
  -F        Also output sequences that have gss.flags != 0.
  -q        Output associated .qual file also.
  -L <shortname>,<shortname>,...
            Only output sequences in libraries with the given shortnames.
  -O common_name,common_name
            Only output sequences in libraries that are from the given organism.
  -E do an SQL EXPLAIN on the generated query, print the results, and exit

  Based loosely on the query-untrimmed-ests.pl from the old EST/Unigene pipeline.
EIEIO

}



usage() unless $ARGV[0] and $ARGV[0] ne "help";

my %opt;
getopts('EFTgql:L:O:',\%opt) or usage();

my ($basename) = @ARGV;

### connect to the local database ###
my $begintime = time();

my $qual_fname = "$basename.qual";
my $seq_fname = "$basename.seq";

open my $seq_out,">$seq_fname"
  or die "Can't open sequence output file \"$seq_fname\" ($!)";

my $qual_out;
if($opt{'q'}) {
  open $qual_out,">$qual_fname"
    or die "Can't open quality output file \"$qual_fname\" ($!)";
}

my $identifier_sql = CXGN::Genomic::GSS->
  external_identifier_sql(qw/ l.shortname
			      c.platenum
			      c.wellrow
			      c.wellcol
			      chr.primer
			      chr.chromat_id
			      g.version
			      /);

$identifier_sql = 'g.gss_id' if $opt{g};

my $seq_sql = do {
  if( $opt{T} ) {
    # untrimmed seq
    "g.seq"
  } else {
    # trimmed seq
    "substring(g.seq from q.hqi_start+1 for q.hqi_length)"
  }
};
my $qual_sql = do {
  if( $opt{q} ) {
    if( $opt{T} ) {
      # untrimmed qual
      ", g.qual"
    } else {
      # trimmed qual
      ", array_to_string( (string_to_array(g.qual,' '))[ q.hqi_start+1 : q.hqi_start+q.hqi_length ] , ' ' )"
    }
  } else {
    ''
  }
};


my $query = <<EOQ;
select $identifier_sql as id,
$seq_sql
$qual_sql
from genomic.library l
join sgn.accession using(accession_id)
join sgn.organism on (sgn.organism.organism_id = sgn.accession.organism_id)
join sgn.common_name cn using(common_name_id)
join genomic.clone c using(library_id)
join genomic.chromat chr using(clone_id)
join genomic.gss g using(chromat_id)
join genomic.qc_report q using(gss_id)
EOQ

my @where;

# -L option
if($opt{L}) {
  #get ids for each of the libs
  my @libids;
  foreach my $ln (split ',',$opt{L}) {
    my @libs = CXGN::Genomic::Library->search(shortname=>$ln);
    @libs == 1
      || die scalar(@libs)." libraries found with shortname '$ln'.\n";
    push @libids,$libs[0]->library_id;
  }

  #    warn 'got libids '.join(' ',keys %libids)."\n";
  @libids or die "Must specify at least one library name with -L.\n";

  push @where, 'c.library_id IN('.join(',',@libids).')';
}

# -F option
unless($opt{F}) {
  push @where,"g.flags = 0";
}

# -l option
if($opt{l}) { 
  die "-l argument must be numeric\n" if $opt{l} =~ /\D/;
  unless($opt{T}) {
    push @where,"q.hqi_length >= $opt{l}";
  } else {
    push @where,"length(g.seq) >= $opt{l}";
  }
}

# -O option
# if($opt{O}) {
#     my $org_sql = join ',', map qq|'$_'|, map lc, split /,/,$opt{O}; #< sql-quote the org names
#     push @where, <<EOQ;
# l.accession_id IN( select a.accession_id
#                    from sgn.accession a
#                    join sgn.organism using(organism_id)
#                    join sgn.common_name cn using(common_name_id)
#                    where lower(cn.common_name) IN($org_sql)
#                  )
# EOQ
# }
if($opt{O}) {
    my $org_sql = join ',', map qq|'$_'|, map lc, split /,/,$opt{O}; #< sql-quote the org names
    push @where, <<EOQ;
lower(cn.common_name) IN($org_sql)
EOQ
}

# # -S option
# if($opt{S}) {
#   my $setstring = join(' | ',(split ',',$opt{S}));
  
#   $gss_query->status_set($setstring);
# }

$query .= "where ".join(' AND ',@where) if @where;

print "Searching for matching sequences, please wait...\n";
my $dbh = CXGN::DB::Connection->new;
#print "$query\n";

if( $opt{E} ) {
    my $a = $dbh->selectall_arrayref("explain $query");
    print map "$_\n", map @$_, @$a;
    $dbh->disconnect(42);
    exit;
}

my $sth = $dbh->prepare($query);
$sth->execute();

print "Dumping sequences ('.' = 100 sequences dumped)...\n";

#write out the sequence and quality data
my $n = 0;
my $empty = 0;
my $empty_qual = 0;

while(my $row = $sth->fetchrow_arrayref) {

    print '.' unless !$n || ($n+$empty) % 100;

    if(!$row->[1]) {
	$empty++;
	next;
    } else {
        $seq_out->print(">$row->[0]\n$row->[1]\n");
    }

    if( $opt{q} ) {
        if( $row->[2] ) {
            $qual_out->print(">$row->[0]\n$row->[2]\n");
        }
        else {
            $empty_qual++;
        }
    }
    $n++;
}

# die "ERROR: Dumped '$n' sequences and skipped '$empty' empties but the query reported finding '".$gss_seqio_orig->result->total_results."' sequences"
#   unless ($n+$empty) == $gss_seqio_orig->result->total_results;

my $runtime = time()-$begintime;
printf("\n%d sequences produced in %d seconds (%.2f seqs/second).\n%d empty sequences were skipped, %d empty qual were skipped.\n",$n,$runtime,$n/($runtime || 1),$empty,$empty_qual);


$dbh->disconnect(42);
