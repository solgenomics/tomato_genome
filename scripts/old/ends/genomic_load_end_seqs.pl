#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;
use Pod::Usage;

use Data::Dumper;

use Bio::PrimarySeq;
use Bio::FeatureIO;

use CXGN::DB::Connection;

use CXGN::BlastDB;
use CXGN::Tools::List qw/ min max any flatten all /;

use CXGN::Genomic::Clone;
use CXGN::Genomic::CloneIdentifiers qw/ parse_clone_ident /;

############## CONFIG AND DEFAULTS #############

our $qual_min           = 20;       # minimum quality value for a 'good' base
our $min_insert_size    = 150;      # minimum size for a high-quality insert to not be flagged as short
our $max_expected_error = 0.035;    # maximum expected error rate tolerable in a sequence
our $min_entropy        = 0.825;    # minimum entropy for a sequence to be considered non-repetitive

############## /CONFIG AND DEFAULTS ############

our %opt;
getopts('q:p:N',\%opt) or pod2usage(1);
$opt{q} ||= $qual_min;
$qual_min =~ /^\d+$/ && $qual_min >= 0
  or die "-q must be a positive integer\n";
$opt{N} ||= 0; #< -N defaults to off

my %primers = parse_cmdline_primer_spec( $opt{p} );
$primers{3} && $primers{5}
  or die "must provide both a 3' and 5' primer with -p\n";

# check for these executables, die if not present
_check_execs_in_path(qw/ vecscreen blastall /);

@ARGV or pod2usage(1);
my ($seq_file,$qual_file) = @ARGV;
-r $seq_file or die "seq file '$seq_file' not found or not readable\n";
!$qual_file || -r $qual_file or die "qual file '$qual_file' not found or not readable\n";

# open our db connection, with AUTOCOMMIT ON
our $dbh = CXGN::DB::Connection->new({ dbargs => {AutoCommit => 0}});

# look up the sequence for the vector in the DB
# my $vector_seq = fetch_vector_seq($opt{V});

my $annots = {};
# hashref as seq_name => { hq => high_quality_feature,
# 			   vec => [vector_feat,...],
# 			   contam => [contam_feat,...],
#                          length => length of seq,
# 		          }

# check that the seqs and quals match up, and also count the length of
# each sequence, and that all the clones are in the database
my $distinct_libs_plates = check_seq_qual($annots,$seq_file,$qual_file,@primers{5,3});
# returns hash as <libname> => { <platenum> => 1 }

# add annotations for vector trimming
annotate_vector_seq($annots,$seq_file);
#hash of seq_name => [vector_feature,...]

# add annotations for quality trimming
annotate_high_quality($annots,$seq_file,$qual_file);
#hash of seq_name => high_quality_feature (Bio::SeqFeature::Generic)

# add annotations for contamination
annotate_contaminants($annots,$seq_file);

# make an efficient where clause to match what clones we're dealing with here
my $clone_whereclause = join ' OR ',
  map {
    my $lib = $_;
    my $min_plate = min( keys %{$distinct_libs_plates->{$lib}} );
    my $max_plate = max( keys %{$distinct_libs_plates->{$lib}} );
    "( l.shortname = '$_' AND c.platenum >= $min_plate AND c.platenum <= $max_plate )"
  } keys %$distinct_libs_plates;
$distinct_libs_plates = undef; #< free up the memory

# load chromatograms
load_chromatograms($seq_file,$clone_whereclause,$primers{5},$primers{3});

# calculate the high-quality insert regions of all the sequences
foreach my $rec (values %$annots) {
  calculate_hqi($rec);
}

# load the gss table
load_gss($annots,$seq_file,$qual_file,$clone_whereclause);

# load the qc_report table
load_qc_report($annots,$seq_file,$clone_whereclause,$primers{5},$primers{3});

# everything looks fine if we've gotten here, commit the load
$dbh->commit;

exit;

############# SUBROUTINES

# return hashref as { seq_name => high_quality_feature }
sub annotate_high_quality {
  my ($annots,$seq_file,$qual_file) = @_;

  warn "annotating high-quality regions...\n";

  my $seq_in  = Bio::SeqIO->new( -file => $seq_file,  -format => 'fasta');

  if( $qual_file ) {
      my $qual_in = Bio::SeqIO->new( -file => $qual_file, -format => 'qual' );
      while ( my $seq = $seq_in->next_seq and my $qual = $qual_in->next_seq ) {

          #get the high-quality region
          #warn "\tfinding quality for ".$qual->primary_id."...\n";
          my $hq = _high_quality_chunk( $qual->qual )
              or next;
          #warn "\tdone.\n";

          #make a feature for it and store it in our global annotation data structure
          $annots->{$seq->primary_id}->{hq} =
              Bio::SeqFeature::Generic->new( -start => $hq->{s} + 1,
                                             -end   => $hq->{e} + 1,
                                             -strand => 0,
                                             -primary => 'high_quality',
                                             -seq_id => $seq->primary_id,
                                           );
      }
      warn "done annotating high-quality regions.\n";
  }
  #### if we have no qual file, just annotate the whole seq as
  #### high-quality
  else {
      while ( my $seq = $seq_in->next_seq ) {
          #make a feature for it and store it in our global annotation data structure
          $annots->{$seq->primary_id}->{hq} =
              Bio::SeqFeature::Generic->new( -start => 1,
                                             -end   => $seq->length,
                                             -strand => 0,
                                             -primary => 'high_quality',
                                             -seq_id => $seq->primary_id,
                                           );
      }
  }
}

#given an arrayref of quality values, find the high-quality chunk in
#it, or return nothing if there is none.  return hashref as
# { s => hq start,
#   e => hq end,
#   length => hq region length,
# }
# note: coordinates are zero-based
sub _high_quality_chunk {
  my ( $qual_scores ) = @_;
  my ( $abs_start,$abs_end ) = ( 0, scalar(@$qual_scores)-1 );
  $abs_start++ while $abs_start < scalar(@$qual_scores) && $qual_scores->[$abs_start] < $qual_min;
  $abs_end--   while $abs_end   > 0                     && $qual_scores->[$abs_end]   < $qual_min;

  return unless $abs_start <= $abs_end;
  #warn "abs_start $abs_start, abs_end $abs_end\n";

  # now quasi-gaussian-filter the qual-scores to smooth out any
  # (very) local quality dips

  my $g_qual_scores = _g_qual_scores_2($qual_scores);
#   my @composite = map {
#     "$qual_scores->[$_]   $g_qual_scores->[$_]"
#   } 0..$#{$qual_scores};
#   die Dumper \@composite;

  #now find the biggest quality-OK chunk in the gaussian-filtered
  #one, then clip it with the absolute min and max
  my $best_chunk = _find_best_qual_chunk($g_qual_scores,$abs_start,$abs_end);

  return $best_chunk if $best_chunk->{length} > 0;
  return;
}

#look at the gaussian-filtered quality scores and find the biggest high-quality chunk
sub _find_best_qual_chunk {
  my ($g_qual_scores,$abs_start,$abs_end) = @_;

  my $last_qual = 0;
  my $best_chunk = { length => 0, s => -1, e => -1};
  my $curr_chunk_begin = $abs_start;
  for (my $idx = $abs_start; $idx <= $abs_end+1; $idx++) {
    #down edge
    if ( !defined($g_qual_scores->[$idx]) || $g_qual_scores->[$idx] < $qual_min && $last_qual > $qual_min  or $idx == $abs_end+1) {
      if ( (my $l = $idx-$curr_chunk_begin) > $best_chunk->{length} ) {
	$best_chunk->{s} = $curr_chunk_begin;
	$best_chunk->{e} = $idx-1;
	$best_chunk->{length} = $l;
      }
    }
    #up edge
    elsif ( $g_qual_scores->[$idx] > $qual_min && $last_qual < $qual_min ) {
      $curr_chunk_begin = $idx;
    }
    $last_qual = $g_qual_scores->[$idx];
  }
  return $best_chunk;
}

#do a quasi-gaussian filter of a quality score array, which smooths out any point defects
sub _g_qual_scores {
  my ($qual_scores) = @_;
  my %filt = ( -2 => 0.15, -1 => 0.2, 0 => 0.3, 1=> 0.2, 2 => 0.15 );
  return [ @{$qual_scores}[0,1],
	   (
	    map {
	      my $idx = $_;
	      my $q = 0;
	      while ( my ($offset,$factor) = each %filt ) {
		$q += $factor * $qual_scores->[$idx+$offset];
	      }
	      $q
	    } (2..scalar(@$qual_scores)-3)
	   ),
	   @{$qual_scores}[-2,-1],
	 ];
}

#NOTE: THIS FUNCTION IS NOT PRETTY, IT IS OPTIMIZED FOR SPEED
sub _g_qual_scores_2 {
  my $qual_scores = shift;
  my @g_filtered_qual = @{$qual_scores}[0,1];
  push @g_filtered_qual,
    map {
	0.15 * $qual_scores->[ $_ - 2 ] +
	0.2  * $qual_scores->[ $_ - 1 ] +
	0.3  * $qual_scores->[ $_     ] +
	0.2  * $qual_scores->[ $_ + 1 ] +
	0.15 * $qual_scores->[ $_ + 2 ];
    } (2..scalar(@$qual_scores)-3);
  push @g_filtered_qual, @{$qual_scores}[-2,-1];
  return \@g_filtered_qual;
}

#add annotations for vector to annotation hash
sub annotate_vector_seq {
  my ($annots,$seq_file) = @_;

  warn "annotating vector sequences...\n";

  my ($univec,$mult) = CXGN::BlastDB->search_like( file_base => '%UniVec' );
  $univec or die "Cannot find UniVec blast DB in ".CXGN::BlastDB->table." table\n";
  $mult and die "Multiple univec databases found in ".CXGN::BlastDB->table." table, I don't know which one to pick";

  sub _vs_parse {
    my $annots = shift;
    my $file = shift;
    my $vec_features = Bio::FeatureIO->new( -file => $file, -format => 'vecscreen_simple' );
    while( my $f = $vec_features->next_feature ) {
      my $seqid = $f->seq_id;
      my $rec = $annots->{$seqid};
      push @{$rec->{vec}}, $f;
      if( $f->end > $rec->{length} ) {
	die "vector is beyond the end of the sequence $seqid, which has length $rec->{length}\n"
            .Dumper($f).Dumper($rec);
      }
    }
  }

  my $vecscreen_save_file = "$seq_file.vecscreen";
  if( $ENV{GENOMIC_LOAD_TESTING} && -f $vecscreen_save_file ) {
    _vs_parse($annots,$vecscreen_save_file);
  } else { 
    my $vecscreen = CXGN::Tools::Run->run( 'vecscreen',
					   -f => 3,
					   -i => $seq_file,
					   -d => $univec->full_file_basename,
					 );
    system 'cp', $vecscreen->out_file, $vecscreen_save_file if $ENV{GENOMIC_LOAD_TESTING};
    eval{ _vs_parse($annots,$vecscreen->out_file) };
    if( $EVAL_ERROR ) {
      system cp => $vecscreen->out_file => $vecscreen_save_file;
      die $EVAL_ERROR."\nvecscreen output saved in $vecscreen_save_file\n";
    }
  }

  warn "done annotating vector sequences.\n";
}

# return hashref as { seq_name => [contam_feat, contam_feat] }
sub annotate_contaminants {
  my ($annots, $seq_file) = @_;

  warn "annotating contaminants...\n";

  my @contam_dbs = map {
    my @o = CXGN::BlastDB->search( file_base => $_ );
    @o == 1 or die "No Blast DB found with file_base $_, maybe you need to update the contamination BDBs used by this script";
    @o
  } qw(
       screening/organelle/tomato_chloroplast
       screening/organelle/ATH1_mitochondria
       screening/lambda-phage/phage-genome
      );
#       screening/rRNA/rRNA-lycopersicon-esculentum

  #  my $seq_in  = Bio::SeqIO->new( -file => $seq_file,  -format => 'fasta');
  foreach my $bdb (@contam_dbs) {
    my $hit_count;

    warn "\tannotating with contaminant '".$bdb->title."'...\n";
    my $contam_save_filename = $bdb->title;
    $contam_save_filename =~ s/\s/_/g;
    $contam_save_filename = "$seq_file.$contam_save_filename.m8";
    if( $ENV{GENOMIC_LOAD_TESTING} && -f $contam_save_filename ) {
      $hit_count = _blast_parse($annots,$contam_save_filename);
    } else { 
      my $blast = CXGN::Tools::Run->run( 'blastall',
					 -i => $seq_file,
					 -e => '1e-10',
					 -p => 'blastn',
					 -m => 8,
					 -d => $bdb->full_file_basename,
				       );
      system 'cp', $blast->out_file, $contam_save_filename if $ENV{GENOMIC_LOAD_TESTING};
      $hit_count = _blast_parse($annots,$blast->out_file);
    }
    warn  "\tdone, found $hit_count contaminated sequences.\n";
  }
  warn "done annotating contaminants\n";
}

#given annots data structure and blast output file, parse it into the data structure
#and return the number of unique sequences that were hit  
sub _blast_parse {
  my $annots = shift;
  my $file = shift;

  open my $blastfh, $file
    or die "$! opening $file";
  my %hits;
  while ( my $line = <$blastfh> ) {
    my ($qname,$hname, $percent_id, $hsp_len, $mismatches, $gapsm,
	$qstart,$qend,$hstart,$hend,$evalue,$bits) = split /\s+/,$line;
    next unless $hsp_len >= 30; #< ignore hits less than 30nt
    $hits{$qname} = 1;
    push @{$annots->{$qname}{contam}},
      Bio::SeqFeature::Generic->new( -start => $qstart,
				     -end   => $qend,
				     -seq_id => $qname,
				     -primary_tag => 'contaminant',
				     -display_name => $hname,
				   );
  }
  return scalar(values %hits);
}

sub _check_execs_in_path {
  foreach (@_) {
    `which $_` or die "$_ executable not found in path, but it's needed for this script to run\n";
  }
}

sub read_class_id {
  my ($classname) = @_;
  our %rc_cache;
  return $rc_cache{$classname} ||= do {
    warn "   looking up read_class_id for $classname\n";
    my $read_class_id = $dbh->selectall_arrayref('select read_class_id from genomic.read_class where lower(class_name) = lower(?)',undef,$classname);
    @$read_class_id or die "no read class found matching '$opt{R}'";
    @$read_class_id > 1 and die scalar(@$read_class_id)." read classes found matching name '$opt{R}' (case insensitive)\n";
    warn "   done.\n";
    $read_class_id->[0][0];
  };
}

# load whatever chromatograms are needed, and make sure they all have clones
sub load_chromatograms {
    my ($seq_file,$clone_whereclause,$primer5,$primer3) = @_;

    warn "loading rows into genomic.chromat table...\n";

    #warn "where $whereclause\n";

    #look up a list of all (relevant) chromats in the DB right now, hash
    #it for efficient lookup, and use this list to decide what we're
    #going to load
    my $chromats_currently_present = $dbh->selectall_hashref(<<EOQ,['shortname','platenum','wellrow','wellcol','primer']);
select
  l.shortname,
  c.clone_id,
  c.platenum,
  c.wellrow,
  c.wellcol,
  chr.primer
from genomic.library l
left join genomic.clone c using(library_id)
left join genomic.chromat chr using(clone_id)
where
  ($clone_whereclause)
EOQ

    #make sure the chromat ID sequence is correctly set
    $dbh->do(<<EOQ);
select setval('genomic.chromat_chromat_id_seq',(select max(chromat_id) from genomic.chromat))
EOQ

    #  die Dumper $chromats_currently_present;
    $dbh->do(<<EOQ);
COPY genomic.chromat
  ( clone_id, primer, direction, read_class_id)
FROM STDIN
EOQ

    open my $idents, "grep '>' $seq_file |"
        or die "$! executing grep on $seq_file\n";
    while (my $ident = <$idents>) {
        chomp $ident;
        ($ident) = $ident =~ /^>(\S+).*$/;

        my $p = parse_clone_ident( $ident, 'bac_end' )
            or die "this should have already been checked!";

        my $dir = $p->{primer} eq $primer5 ? 5 :
            $p->{primer} eq $primer3 ? 3 :
                die "don't know whether primer '$p->{primer}' is 5' or 3'";

        my $read_class_id = read_class_id("$p->{clonetype} ends");

        my $clone_rec = $chromats_currently_present->{$p->{lib}}{$p->{plate}}{$p->{row}}{$p->{col}}
            or die "no clone found for $ident, are all clones loaded?\n".Dumper($p);
        if (!$opt{N} && ref $clone_rec->{$p->{primer}} ) {
            warn "chromatogram $ident already loaded.\n";
        } else {
            my ($first_primer) = keys %$clone_rec;
            $dbh->pg_putline( join("\t",
                                   ($clone_rec->{$first_primer}{clone_id} or die('invalid clone rec: '.Dumper $clone_rec)),
                                   $p->{primer},
                                   $dir,
                                   $read_class_id,
                                  )
                              ."\n"
                            );
        }
        #    print( join("\t",($clone_rec->{''}{clone_id} or die( Dumper $clone_rec)),uc($p->{primer}),$dir,$read_class_id)."\n" );
    }
    $dbh->pg_endcopy;

    warn "done loading genomic.chromat\n";
}

# load necessary rows in gss table
sub load_gss {
  my ($annots,$seq_file,$qual_file,$clone_whereclause) = @_;

  warn "loading rows into genomic.gss table...\n";

  #query which relevant gss's already exist in the database
  #going to load
  my $gss_currently_present = $dbh->selectall_hashref(<<EOQ,['shortname','platenum','wellrow','wellcol','primer','version']);
select
  l.shortname,
  c.clone_id,
  c.platenum,
  c.wellrow,
  c.wellcol,
  chr.chromat_id,
  chr.primer,
  g.version,
  g.gss_id,
  g.seq,
  g.qual
from genomic.library l
left join genomic.clone c using(library_id)
left join genomic.chromat chr using(clone_id)
left join genomic.gss g using(chromat_id)
where
  ($clone_whereclause)
order by chr.date
EOQ

  my $seq_in  = Bio::SeqIO->new( -file => $seq_file,  -format => 'fasta');
  my $qual_in = $qual_file ? Bio::SeqIO->new( -file => $qual_file, -format => 'qual' ) : undef;

#make sure the gss_id seq is correctly set
  $dbh->do(<<EOQ);
select setval('genomic.gss_gss_id_seq',(select max(gss_id) from genomic.gss))
EOQ
  my $copy = $dbh->do(<<EOQ);
 COPY genomic.gss
   ( chromat_id, version, seq, qual, status, flags )
 FROM STDIN
EOQ

  while( my $seq = $seq_in->next_seq ) {

    my $ident = $seq->primary_id;

    my $qual;
    my $qual_text;
    if( $qual_file ) {
        $qual = $qual_in->next_seq or die "what??? end of qual file!\n";
        $qual_text = join ' ',@{$qual->qual};

        #make sure we have the right qual
        $ident eq $qual->primary_id
            or die "quals and seqs not in same order!\n";
    }

    my $p = parse_clone_ident( $ident, 'bac_end' )
      or die "this should have already been checked!";

    my $chromat_rec = $gss_currently_present->{$p->{lib}}{$p->{plate}}{$p->{row}}{$p->{col}}{$p->{primer}}
      or die "no chromat found for $ident, and it should be there by now.\nclone end id parse dump: ".Dumper($p);

    #is there an existing version
    my $max_existing_version = max( keys %$chromat_rec );
    my $new_version = ($max_existing_version||0) + 1;

    my $seqannots = $annots->{$ident} || {}; #< the annotations entry for this sequence

    #calculate entropy and error BEFORE we decide whether to load it, because we might
    #need the entropy and error for the qc_report table
    if( $seqannots->{hqi} ) { #<if we have an HQ insert
      $seqannots->{hqi_entropy} = find_entropy($seq->subseq($seqannots->{hqi}->location));
      $seqannots->{hqi_error}   = find_expected_error( $seqannots->{hqi}->start-1,
						       $seqannots->{hqi}->length,
						       ($qual ? $qual->qual : undef),
						     );
    }

    #if there are existing GSS's for this chromat, check whether this
    #one is really different, and warn and skip it if not
    my $newest_existing_gss = $chromat_rec->{$max_existing_version};
    if( $newest_existing_gss->{gss_id}
	&& $newest_existing_gss->{seq} eq $seq->seq
	&& ( !defined $newest_existing_gss->{qual} && !defined $qual_text
             || $newest_existing_gss->{qual} eq $qual_text
           )
      ) {
      warn "skipping load for $ident already exists with gss id $newest_existing_gss->{gss_id}\n";
      next;
    }

    my $quality_flags = 0;

    if( $seqannots->{hqi} ) { #<if we have an HQ insert

      #check for short insert size
      $quality_flags |= 0x4 unless $seqannots->{hqi}->length >= $min_insert_size;

      #check for high expected error
      $quality_flags |= 0x8 unless
	$max_expected_error > $seqannots->{hqi_error};

      #check for low entropy
      $quality_flags |= 0x10 unless $min_entropy < $seqannots->{hqi_entropy};

    } else {#< if we don't have an HQ insert
      $quality_flags |= 0x4; #< set short insert since there appears to be no insert
    }


    #check for any contamination features regardless of whether they
    #overlap with the HQ region
    $quality_flags |= 0x20 if $seqannots->{contam} && @{$seqannots->{contam}};

    $dbh->pg_putline( join("\t",
                           ($chromat_rec->{$max_existing_version}{chromat_id} or die('invalid chromat rec: '.Dumper $chromat_rec)),
			   $new_version,
			   $seq->seq,
			   $qual_text || '\N',
			   #status_flags
			   ( 0x40 # chimera_unk
			     | 0x80 # repeats_unk
			   ),
			   $quality_flags,
			  )
		      ."\n"
		    );
  }
  $dbh->pg_endcopy;
  warn "done loading genomic.gss\n";
}


#now load any necessary rows into qc_report
sub load_qc_report {
  my ($annots,$seq_file,$clone_whereclause,$primer5,$primer3) = @_;

  warn "loading rows into genomic.qc_report table...\n";

  #query which relevant gss's already exist in the database
  #going to load
  my $qcr_currently_present = $dbh->selectall_hashref(<<EOQ,['shortname','platenum','wellrow','wellcol','primer','version']);
select
  l.shortname,
  c.platenum,
  c.wellrow,
  c.wellcol,
  chr.primer,
  g.version,
  g.gss_id,
  q.qc_report_id
from genomic.library l
left join genomic.clone c using(library_id)
left join genomic.chromat chr using(clone_id)
left join genomic.gss g using(chromat_id)
left join genomic.qc_report q using(gss_id)
where
  ($clone_whereclause)
order by chr.date
EOQ

#make sure the qc_report_id seq is correctly set
  $dbh->do(<<EOQ);
select setval('genomic.qc_report_qc_report_id_seq',(select max(qc_report_id) from genomic.qc_report))
EOQ
  my $copy = $dbh->do(<<EOQ);
 COPY genomic.qc_report
   ( gss_id, vs_status, qstart, qend, hqi_start, hqi_length, entropy, expected_error, qual_trim_threshold )
 FROM STDIN
EOQ

  open my $idents, "grep '>' $seq_file |"
    or die "$! executing grep on $seq_file\n";
  while(my $ident = <$idents>) {
    chomp $ident;
    ($ident) = $ident =~ /^>(\S+)/;

    my $p = parse_clone_ident( $ident, 'bac_end' )
      or die "this should have already been checked!";

    my $chromat_rec = $qcr_currently_present->{$p->{lib}}{$p->{plate}}{$p->{row}}{$p->{col}}{$p->{primer}}
      or die "no chromat found for $ident, and it should be there by now.\nclone end id parse dump: ".Dumper($p);

    #is there an existing GSS version
    my $max_existing_version = max( keys %$chromat_rec );
    $max_existing_version or die "no GSS row found, and there should be one by now.\n".Dumper($chromat_rec);;

    my $gss_rec = $chromat_rec->{$max_existing_version};
    $gss_rec->{gss_id} or die "sanity check failed";

    #if there are existing QCRs's for this GSS, check whether this
    #one is really different, and warn and skip it if not
    if( my $existing_qcr = $gss_rec->{qc_report_id} ) {
      warn "skipping existing qcreport $gss_rec->{qc_report_id} for gss $gss_rec->{gss_id}\n";
      next;
    }

    my $seqannots = $annots->{$ident} || {}; #< the annotations entry for this sequence

    my $vs_status = do {
      my $got_vec_flag      = $seqannots->{vec} && @{$seqannots->{vec}};
      my $no_insert_flag    = !$seqannots->{hqi};
      my $short_insert_flag = $seqannots->{hqi} && $seqannots->{hqi}->length <= $min_insert_size;

      my $dir = $p->{primer} eq $primer5 ? 5 :
	        $p->{primer} eq $primer3 ? 3 :
		  die "don't know whether primer '$p->{primer}' is 5' or 3'";

      $no_insert_flag      ? "noinsert$dir" :
	$short_insert_flag ? "short$dir"    :
	  $got_vec_flag    ? "good$dir"     :
	                     'novec';
    };
    $vs_status = CXGN::Genomic::QCReport::qc_num($vs_status);
    defined $vs_status or die("cannot find status number for '$vs_status'"),

    my ($qstart,$qend) = (0,0);
    ($qstart,$qend) = ($seqannots->{hq}->start-1, $seqannots->{hq}->end - 1) if $seqannots->{hq};

    # NOTE: we are not loading istart and iend, they're not really important

    my ($hqistart,$hqilength) = (0,0);
    ($hqistart,$hqilength) = ($seqannots->{hqi}->start-1, $seqannots->{hqi}->length) if $seqannots->{hqi};

    $dbh->pg_putline( join( "\t",
			    $gss_rec->{gss_id},
			    $vs_status,
			    $qstart,$qend,
			    $hqistart,$hqilength,
			    $seqannots->{hqi_entropy} || '\N',
			    $seqannots->{hqi_error} || '\N',
			    $qual_min
			  )
		      ."\n"
		    );
  }
  $dbh->pg_endcopy;
  warn "done loading genomic.qc_report\n";
}

# return a Bio::SeqFeature::Generic of the region of the sequence
# that is both high-quality and does not have a match to vector
# sequence
sub calculate_hqi {
  my ($seqannot) = @_;

  return unless $seqannot->{hq};

  my $hq = $seqannot->{hq}; #< the high-quality region in this seq
  if( $hq->end > $seqannot->{length} ) {
    die "invalid hq for sequence with length $seqannot->{length}: ".Dumper($hq);
  }

  my $seqid = $hq->seq_id;
  my $hqi = Bio::SeqFeature::Generic->new( -start  => $hq->start,
					   -end    => $hq->end,
					   -seq_id => $hq->seq_id,
					 );

  #now go through the cut down the HQI with each vector match
  foreach my $vec (@{$seqannot->{vec}}) {

    my ($vs,$ve) = ($vec->start,$vec->end);
    my ($hs,$he) = ($hqi->start,$hqi->end);

    #if our hqi isn't valid, it must be all gone
    return if $hs < 1 || $he > $seqannot->{length} || $hs > $he;

    $vs && $ve && $hs && $he
      or die "not all there: $vs, $ve, $hs, $he : ".Dumper($seqannot);

    if( $ve > $seqannot->{length} ) {
      die "vector region $vs,$ve is beyond the end of the sequence $seqid, which has length $seqannot->{length}\n".Dumper($seqannot);
    }

    #in this case, vec completely overlays remaining hqi, so just return no hqi
    return if $vs <= $hs && $ve >= $he;

    #in this case, vec lies outside of remaining hqi, so does not affect it
    next if $vs > $he || $ve < $hs;

    # in this case, vec is completely internal to hqi, so remaining hqi
    # is the biggest piece left
    if( $vs > $hs && $ve < $he ) {
      my $left_length = $vs - $hs;
      my $right_length = $he - $ve;
      if( $right_length > $left_length ) { #< take the right piece
	$hqi->start( $ve+1 );
      } else { #< take the left piece if they are equal or the left is bigger
	$hqi->end( $vs-1 );
      }
      next;
    }

    #in this case, vec contains left end
    if( $ve >= $hs ) {
      $hqi->start( $ve+1 );
      next;
    }

    #in this case, vec contains right end
    if( $vs <= $he ) {
      $hqi->end( $vs-1 );
      next;
    }
  }

  #return nothing if we've covered the whole hq region with vector
  return unless $hqi->start <= $hqi->end;
  return $seqannot->{hqi} = $hqi;
}

#given a vector name, fetch it from the DB or die if it's not there
sub fetch_vector_seq {
  my ($vecname) = @_;

  my $v = $dbh->selectall_arrayref(<<EOQ,undef,$vecname);
select name,seq
from sgn.cloning_vector
where lower(name) = lower(?)
EOQ

  @$v or die "no vector found in sgn.cloning_vector table named '$vecname'\n";
  @$v > 1
    and die "multiple matching vectors found:\n",
      map "   $_->[0]\n", @$v;

  return Bio::PrimarySeq->new( -primary_id => $v->[0][0],
			       -seq => $v->[0][1]
			     );
}

# read in and match up the seq and qual files, making sure the identifiers are parsable,
# and read the lengths of all the sequences
# die on error
sub check_seq_qual {
  my ($annots, $seq_file, $qual_file, $primer5, $primer3) = @_;

  warn "validating sequences and qual input data...\n";

  our $seq_error_flag = 0;
  sub _seq_error($) {
    warn @_,"\n";
    $seq_error_flag = 1;
  }

  die "Errors found in seq and qual files.  Aborting.\n" if $seq_error_flag;

  my %distinct;

  my $seq_in = Bio::SeqIO->new( -file => $seq_file, -format => 'fasta');
  my $qual_in;
  if( $qual_file ) {
      $qual_in = Bio::SeqIO->new( -file => $qual_file, -format => 'qual' );
  }

  while ( my $seq = $seq_in->next_seq ) {

    my $seqid = $seq->primary_id;


    my $qual;
    if( $qual_file ) {
        $qual = $qual_in->next_seq
            or die "Qual and seq file do not have the same number of sequences in them.  More seqs than quals.\n";
        $qual->primary_id eq $seqid
            or die "Qual and seq files do not have the same sequences in the same order.  Maybe you should sort them by seq name and try again?\n";
    }

    my $l = $seq->length;
    if( defined $annots->{$seqid}{length} ) {
      _seq_error("$seqid appears multiple times in seq and qual files!  Loading multiple reads from the same clone end requires multiple runs of this script.\n");
    } else {
      $annots->{$seqid}{length} = $l;
    }
    unless(!$qual || $qual->length == $l) {
      _seq_error("Seq and qual for ".$seq->primary_id." are not the same length.  (seq: ".$seq->length.", qual: ".$qual->length.")");
    }

    if( $seq->seq =~ /([^ACGTURYKMSWBDHVN]+)/i ) {
      _seq_error("Invalid character(s) found in seq '".$seq->primary_id."': '$MATCH'");
    }

    my $p = parse_clone_ident( $seq->primary_id, 'bac_end' );

    unless($p) {
      _seq_error("Malformed bac end identifier '".$seq->id."'");
    } else {
      $p->{chromat_id} == 0
         or  _seq_error($seq->id." should end in _0, not _$p->{chromat_id}");

      #look up and cache the read_class_id for this type of clone
      my $read_class_id = read_class_id("$p->{clonetype} ends");

      #look up this clone in the database, make sure it's in there
      my $clone = CXGN::Genomic::Clone->retrieve_from_parsed_name($p)
	or _seq_error("no clone seems to be loaded for ".$seq->id);
    }

    unless($p->{primer} eq $primer5 || $p->{primer} eq $primer3) {
      _seq_error("Unknown primer '$p->{primer}' in seq ident '".$seq->primary_id."'.  Do you need to add $p->{primer} to the -p option passed to this script?");
    }

    $distinct{$p->{lib}}{$p->{plate}} = 1;
  }
  if( $qual_file && $qual_in->next_seq) {
    die "Qual and seq files do not have the same number of sequences in them.  More quals than seqs.\n";
  }

  if($seq_error_flag) {
    die "Errors found in sequence and/or qual files.  Aborting.\n";
  }

  warn "done checking seq and qual data.\n";
  return \%distinct;
}

#find the entropy value of a sequence string
sub find_entropy {
    my ($seq) = @_;
    my $log2 = log(2);

    return 0.0 if length($seq) < 32;
    my @seq = split//,$seq;

    my %count = ();
    my $basis = 0;
    for(my $i=0;$i<@seq-1;$i+=2) {
	$count{$seq[$i] . $seq[$i+1]}++;
	$basis++;
    }

    my $e = 0.0;
    foreach ( values %count ) {
	my $p = $_ / $basis;
	$e += $p * log($p)/$log2;
    }
    $e /= -4.0;

    return $e;
}

# find the expected error of the hqi region of an arrayref of qual
# values
sub find_expected_error {
    my ($qs, $ql, $qualref) = @_;

    return 0 unless $qualref;

    return 1.0 if $ql < 2;

    my $ee = 0.0;
    for( my $i=$qs-1; $i<$qs+$ql; $i++ ) {
	$ee += 10**($qualref->[$i]/-10.0);
    }

    return $ee/($ql);
}

# arg: spec string passed as -p arg
# return: list of ( 5 => primer name, 3 => primer name )
sub parse_cmdline_primer_spec {
    my $spec = shift || '';

    my @f = split /,/,$spec;
    @f = map [split /=/,$_], @f;
    pod2usage('invalid primer spec') unless
        @f == 2
            && all( map @$_==2,@f )
                && all( map $_->[0] =~ /^\d+$/, @f );

    return flatten @f;
}

__END__

=head1 NAME

genomic_load_end_seqs.pl - load clone-end sequences

=head1 SYNOPSIS

  genomic_load_end_seqs.pl -p <primer spec> seq_file qual_file

  Load clone-end sequences into the Genomic database.
  Seq file and qual file must have the same identifiers,
  and be in the same order.

  Identifiers must be formatted as

      <SOL-style clone_name>_<primer name>_<num>

      Note that <num> in the above it not used by this loading script,
      the chromat_id will be chosen automatically.  You can just make
      it 0.  Example: C<SL_EcoRI0123A12_SP6_0>

  Options:

    -q <num>
       Quality minimum for determining high-quality region.
       Default: $qual_min

    -p <primer spec>
       Which primer name corresponds to a 5' read, and which
       a 3' read.  Formatted as 5=primer_name,3=primer_name

    -N treat sequences as new reads (chromatograms) which have not
       been loaded before.  Otherwise, treats input seqs as new
       (superseding) basecalls of the same chromatograms if a
       chromatogram already exists for that clone and primer.

  Example:
     $FindBin::Script -p 5=T7,3=SP6 myseqs.seq myseqs.qual

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
