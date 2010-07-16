#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;  $SIG{__DIE__} = \&Carp::confess;
use Data::Dumper;
use FindBin;
use Getopt::Std;
use POSIX;

use File::Basename qw/ basename /;
use File::chdir;
use File::Copy;
use File::Spec;
use File::Path qw/ rmtree mkpath /;
use File::Temp qw/ tempfile tempdir /;
use List::Util qw/ sum min max /;
use List::MoreUtils qw/ uniq any /;
use YAML::Any;

use Bio::Index::Fasta;

use CXGN::Cluster;
use CXGN::Cview::MapFactory;
use CXGN::DB::Connection;

use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident assemble_clone_ident/;

use CXGN::Publish qw/published_as publish/;

use CXGN::Tools::Run;
use CXGN::Tools::Script qw/lock_script unlock_script/;

use CXGN::TomatoGenome::BACPublish qw/aggregate_filename agp_file/;
use CXGN::TomatoGenome::Config;

########### CONFIGURATION/DEFAULTS ################

my $cfg = CXGN::TomatoGenome::Config->load_locked;
my $country_uploads_path = $cfg->{'country_uploads_path'};
my $publish_path = File::Spec->catdir(@{$cfg}{  'ftpsite_root',  'bac_publish_subdir' });
my $agp_path     = File::Spec->catdir(@{$cfg}{  'ftpsite_root',  'agp_publish_subdir' });

#mummer params
my $mummer_min_overlap = 1500;

my $cview_physical_map_version = 'p9';

###################################################

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options]

  Takes the BAC sequences on each chromosome, tries to assemble them
  with mummer, and constructs an AGP file from the results.

  Options:

    -C run mummer jobs on the cluster, not locally.
       WARNING: mummer currently just segfaults when run on cluster
       nodes

    -p <dir>
       set BAC publishing dir to read sequences from.
       Default: $publish_path

    -a <dir>
       set directory for outputting AGP files
       Default: $agp_path

    -l <overlap>
       set the minimum overlap to use with mummer
       Default: $mummer_min_overlap

    -m <map id>
       the Cview map ID to use for establishing the positions of BACs
       Default: $cview_physical_map_version

    -c <list>
       list of chromosome numbers to process.  Default 0..12

    -A if passed, will also publish a complete_assembly.tar.gz file
       containing all phrap output, .ace files, contig membership
       listings, phrap consensus sequences, and generated AGP filesa

EOU
}
sub HELP_MESSAGE {usage()}

my @command_line_args = @ARGV; #< save argv for the metadata file
our %opt;
getopts('Cp:a:m:c:A',\%opt) or usage();
@ARGV and usage(); #< there should be no non-option arguments

#get our publishing path
$publish_path = $opt{p} if defined $opt{p};
-r $publish_path or die "publish path $publish_path not found or not readable\n";

#get our AGP path
$agp_path = $opt{a} if defined $opt{a};
-w $agp_path or die "agp path $agp_path not found or not writable\n";

#get our mummer min overlap
$mummer_min_overlap = $opt{l} if defined $opt{l};
$mummer_min_overlap =~ /^\d+$/ && $mummer_min_overlap > 0
  or die "invalid -l min overlap '$mummer_min_overlap', must be a positive integer\n";
my @mummer_options = (
    '-mum',
    '-b',
    '-n',
    -l => $mummer_min_overlap,
);


#get our chromosome numbers
my @chromosome_nums = (0..12);
@chromosome_nums = eval "($opt{c})" if defined $opt{c};
die "invalid chromosome numbers expression\n" if $EVAL_ERROR;

#get our cview physical map version
$cview_physical_map_version = $opt{m} if $opt{m};

$opt{A} &&= File::Spec->catdir( tempdir( CLEANUP => 1 ), 'bac_assembly' );

# init the map data model
my $dbh = CXGN::DB::Connection->new({ config => $cfg });
my $mf =  CXGN::Cview::MapFactory->new( $dbh, $cfg );
my $physical_map = $mf->create( { map_version_id => $cview_physical_map_version } )
    or die "failed to create physical map object for map version $cview_physical_map_version\n";

lock_script() or die "Do not run more than one $FindBin::Script at the same time.\n";

my %chrdata; #< big hash of all the data about a given chromosome

# find all the BAC sequences to use, and make per-chromosome sequence
# files to use as input to mummer
my $bacs_file = published_as( aggregate_filename("all_seqs",$publish_path) )
    or die "cannot find bacs file at ".aggregate_filename( 'all_seqs', $publish_path );
my $bacs_seqio = Bio::SeqIO->new( -file => $bacs_file->{fullpath}, -format => 'fasta' );

for (@chromosome_nums) {
    my ( $fh,$file ) = tempfile( UNLINK => 1 );
    $chrdata{$_} = { seqfile => $file, seqfile_fh => $fh,  chrnum => $_ };
}

while( my $seq = $bacs_seqio->next_seq ) {
    my $p = parse_clone_ident( $seq->id, 'versioned_bac_seq' )
        or die "could not parse seq id ".$seq->id;
    $p->{chr} = 0 if $p->{chr} eq 'unmapped';
    if( my $chr_rec = $chrdata{ $p->{chr} } ) {
        $chr_rec->{seqfile_fh}->print( ">", $seq->id, "\n", $seq->seq, "\n" );
    } else {
        #warn "skipping seq ".$seq->id." (chr $p->{chr})\n";
    }
}
$_->{seqfile_fh}->close for values %chrdata;

# also, index all the bac seqs we are dealing with
my $seqs_index = Bio::Index::Fasta->new(
    -filename => do { my (undef,$tf) = tempfile(UNLINK => 1); $tf },
    -write_flag => 1
   );
$seqs_index->make_index( $bacs_file->{fullpath} );

# fetch all the known map positions of BACs.  to include a BAC cluster
# in our AGP, we have to know its location on the physical map
fetch_mapping_data( \%chrdata );

#dispatch mummer jobs for each of them
my $runfunc = do {
  if($opt{C}) {
    CXGN::Tools::Run->temp_base('/data/shared/tmp');
    'run_cluster'
  } else {
    'run'
  }
};

foreach my $chr_num ( sort {$a <=> $b} keys %chrdata ) {
    my $chr_rec = $chrdata{$chr_num};
    local $CWD = File::Spec->tmpdir;

    $chr_rec->{job} =
        CXGN::Tools::Run->$runfunc(

            'mummer',

            @mummer_options,

            -F => $chr_rec->{seqfile},

            $chr_rec->{seqfile},

            { die_on_destroy => 1,
              working_dir => File::Spec->tmpdir,
              on_completion => sub {
                  my $job = shift;
                  my $cluster_set  = mummer_to_clusterset( $job->out_file );
                  my @clusters      = map_clusters( $cluster_set, $chr_rec->{mapped_bacs} );
                  #warn "got ".scalar(@clusters)." mapped clusters\n";
                  my $agp_results = generate_agp_file( $chr_rec->{chrnum}, \@clusters );
                  $chr_rec->{agp_file} = $agp_results->{file};
                  $chr_rec->{metadata} = $agp_results->{metadata};
              },
          }
           )
              or die "failed to run mummer on $chr_rec->{seqfile}!";

}

#when the script ends, clean up all the cluster job tempfiles
END { $_->{job} && $_->{job}->cleanup foreach values %chrdata }

my @agp_publish_cmds;

if( $opt{A} ) {

    # write the assembled-contig sequence sets
    make_assembly_dir_contig_files( $opt{A}, $seqs_index );

    # also make a subdir in the assembly dir with the AGP files we generated
    copy_agp_files_to_assembly_dir( $opt{A}, \%chrdata );

    my %build_metadata = map {
        my $v = $chrdata{$_}{'metadata'};
        "chromosome $_ AGP" => $v
    } keys %chrdata;

    write_assembly_metadata( $opt{A},
        'BAC sequence set'       => $bacs_file->{fullpath},
        'chromosomes included in this assembly'
                                 => join( ', ', @chromosome_nums ),
        'mummer options'         => join( ' ', @mummer_options ),
        'phrap options'          => join( ' ', CXGN::Cluster::Precluster->phrap_options ),
        'chromosome AGP builds'  => \%build_metadata,
    );

    my $assembly_tar = File::Temp->new;
    make_complete_assembly_tar( $opt{A}, $assembly_tar );
    push @agp_publish_cmds,
        ['cp', $assembly_tar, File::Spec->catfile( $agp_path, 'bac_assembly.tar.gz' ) ];
}


# we will replace the AGP files only for the chromosomes that don't
# have a manually uploaded one
push @agp_publish_cmds,
    map {
        -f $chrdata{$_}{agp_file} or warn "$chrdata{$_}{agp_file} is not where I left it!\n";
        [ 'cp', $chrdata{$_}{agp_file}, agp_file($_,1,$agp_path) ]
    }
    # list of chr nums that have either no agp, or an autogenerated one
    grep {
        if( my $agp = agp_file($_,0,$agp_path) ) {
            #is the file there autogenerated?  if so, regenerate, if not, leave it alone
            `grep autogenerated-by $agp` ? 'generate' : ''
        } else {
            #< no file at all
            'generate'
        }
    } @chromosome_nums;


publish(@agp_publish_cmds);

unlock_script();

$dbh->disconnect(42);

# system 'cat', map $_->{agp_file},values %chrdata;
# system 'validate_agp.pl', $_->{agp_file} foreach values %chrdata;

################ SUBROUTINES ##########################################################33

sub generate_agp_file {
  my ($chrnum, $contigs) = @_;
#   $Data::Dumper::Maxdepth = 5;
#   warn Dumper $contigs;

  my ($agp_fh,$agp_file) = tempfile(UNLINK => 1);
  print $agp_fh "#Draft Tomato chromosome $chrnum\n";
  print $agp_fh "#autogenerated-by SGN mummer/phrap AGP generator ($FindBin::Script)\n";


  if( $opt{A} ) {
      mkpath( $opt{A} );
      rmtree( $_ ) for glob "$opt{A}/chr${chrnum}_*";
  }

  my %metadata; #< keep some statistics and metadata about this AGP generation

  # we assume that 1cM is on average about 750kB of sequence (Wing, Zhang, and Tanksley, 1994).
  #
  my $correspondence = 100_000;

  my $line_count = 0;
  my $previous_global_end = 0; #< the global end of the previous line
  my $printline = sub(@) {
    my ($s,$e,@other) = @_;
    push @other,'' unless @other == 5;
    ++$line_count;
    print $agp_fh join("\t", 'S.lycopersicum-chr'.$chrnum, $s, $e, $line_count, @other )."\n";
  };


  my $members_fh;
  if( $opt{A} ) {
      # dump the cluster's member list also
      my $members_file = _members_filename( $opt{A} );
      open $members_fh, '>>', $members_file or die "$! writing $members_file";
  }

  my $printed_unmapped_divider = 0; #<flag of whether we have printed the comment about unmapped sequences already
  my $contig_number = 0;
  for( my $precluster_number = 1; $precluster_number <= @$contigs; $precluster_number++ ) {
    my $mapped_contig = $contigs->[$precluster_number-1];
    my $offset = $mapped_contig->{offset};
    if( ! $printed_unmapped_divider && $offset > 500 ) {
      print $agp_fh "# END OF DRAFT CHROMOSOME BUILD.  THE LINES BELOW CONTAIN SEQUENCES THAT COULD NOT BE LOCATED ON THE PHYSICAL MAP\n";
      $printed_unmapped_divider = 1;
    }
    my $cluster = $mapped_contig->{cluster};
    #warn "new cluster: ".join(' ',$cluster->get_members)."\n";
    my $base = sprintf('%0.0f',$offset * $correspondence + 1);
    if( $base <= $previous_global_end+50_001 ) {
      #this means the BACs are mapped right on top of eachother, but they don't have any sequence overlap, so introduce a 'contig no' gap
      $base = $previous_global_end+50_001;
    }
#     if ($mapped_clusters{$cluster_positions{$offset}->get_unique_key()} > 1) { 
#       # determine order of bacs in the contig... TBD
#     }

    #put in an inter-contig gap if necessary
    if($base > $previous_global_end+1) {
      my $gap_start = $previous_global_end+1;
      my $gap_length = $base-$gap_start;
      my $gap_end = $gap_start+$gap_length-1;
      $printline->($gap_start,$gap_end,'N',$gap_length,'contig','no');
      $previous_global_end = $gap_end;
    } elsif( $base < $previous_global_end+1 ) {
      #or adjust the base up if necessary
      $base = $previous_global_end+1;
    }


    # if we are producing an assembly results dir, do stuff for that.
    if( $opt{A} ) {
        my $dir = _precluster_dir( $opt{A}, "chr$chrnum", $precluster_number );
        $cluster->set_assembly_dir( $dir );
    }

    my $contigs_in_precluster = 0;
    foreach my $contig (  $cluster->get_consensus_base_segments( $seqs_index, min_segment_size => 2000 ) ) {
      my $previous_contig_end = 0; #< the contig end of the previous line
      ++$contig_number;

      #put in a 20kb clone gap if necessary, if we have members in a
      #cluster that did not assemble into one contig
      if( ++$contigs_in_precluster > 1 ) {
	my $gap_start = max($base,$previous_global_end+1);
	my $gap_end = $gap_start + 50_000 - 1;
	$printline->( $gap_start, $gap_end, 'N', 50_000, 'clone', 'yes' );
	$base = $gap_end + 1;
	$previous_global_end = $gap_end;
      }

      $members_fh->print( _precluster_name( "chr$chrnum", $precluster_number), "\t", _contig_name( "chr$chrnum", $contig_number ) ) if $members_fh;

      foreach my $member ( @$contig ) {
          my ( $member_contig_start, $member_contig_end, $name, $member_local_start, $member_local_end, $member_reverse ) = @$member;

          $members_fh->print( "\t$name" ) if $members_fh;

          my $seq = $seqs_index->fetch($name)
              or die "could not get seq for '$name'";
          my $seqlength = $seq->length
              or die "no seq length returned for '$name' seq";

          # for AGP, the member coordinates are relative to the
          # *uncomplimented* sequence.  so, for complemented component
          # seqs, switch those around.
          if( $member_reverse ) {
              ( $member_local_start, $member_local_end ) =
                   ( $seqlength - $member_local_end   + 1, $seqlength - $member_local_start + 1 );
          }

          my $member_global_start = $base + $member_contig_start - 1;
          my $member_global_end   = $base + $member_contig_end   - 1;

          # sometimes phrap shortens or lengthens runs of repetitive
          # nucleotides in order to build better consensus sequences.
          # this is usually only a few nt per 100kb.  these edits are
          # very difficult to represent in an AGP, so for now, just
          # ignore edits, fixing up coordinates if they go beyond the
          # end of the sequence.
          #
          # this inability to represent the phrap sequence edits is
          # going to lower the quality of some contigs
          if( $member_local_end > $seqlength ) {
              my $difference = $member_local_end - $seqlength;
              push @{$metadata{'reads with edits introduced by phrap'}}, { $name => "$difference additional bp" };
              warn "WARNING: artificially shortening consensus segment $name ( $member_local_start, $member_local_end ) by $difference bases to cope with phrap sequence edit. This will make a slight error in chr $chrnum contig $contig_number (precluster $precluster_number).\n";
              $member_local_end -= $difference;
              unless( $member_reverse ) {
                  $member_global_end -= $difference;
              }
          }
          elsif( $member_local_start < 1 ) {
              my $difference = 1 - $member_local_start;
              push @{$metadata{'reads with edits introduced by phrap'}}, { $name => "$difference additional bp" };
              warn "WARNING: artificially shortening consensus segment $name ( $member_local_start, $member_local_end ) by $difference bases to cope with phrap sequence edit. This will make a slight error in chr $chrnum contig $contig_number (precluster $precluster_number).\n";
              $member_local_start += $difference;
              if( $member_reverse ) {
                  $member_global_end  -= $difference;
              }
          }

          $printline->( $member_global_start, $member_global_end, 'F',
                        $name, $member_local_start, $member_local_end,
                        $member_reverse ? '-' : '+',
                       );
          $previous_global_end = $member_global_end;
          $previous_contig_end = $member_contig_end;
      }
      $members_fh->print( "\n" ) if $members_fh;
    }
  }

  close $agp_fh; #< gotta close it to flush the buffers

  return { file => $agp_file, metadata => \%metadata };
}

# fetches the mapping data from Cview
sub fetch_mapping_data {
    my $chrdata = shift;

    while ( my ($chrnum,$chr_rec) = each %$chrdata ) {

        next unless $chrnum > 0;

        #  warn "fetching chr for $chrnum\n";
        my $chr_map = $physical_map->get_chromosome($chrnum)
            or die "failed to fetch physical map chromosome $chrnum for map version $cview_physical_map_version\n";

        # make a hash of name => offset of mapped BACs on this chromosome
        my %b =
            map  { $_->get_marker_name => $_->get_offset }
                grep $_->isa('CXGN::Cview::Marker::Physical'),
                    $chr_map->get_markers;

        $chr_rec->{mapped_bacs} = \%b;
    }
}

#$cluster_set is a CXGN::Cluster::ClusterSet
#$map_positions is a hashref of { bac_name => offset }
#returns a list of { cluster => obj, position => num }
#sorted ascending by position
sub map_clusters {
  my ($cluster_set,$map_positions) = @_;

  return
    sort {$a->{offset} <=> $b->{offset}}
    map {
      my $cluster = $_;
      #      warn "in cluster with ".join(' ',$cluster->get_members)."\n";
      #if the cluster has a map position, then we'll use it in the AGP
      if (my @map_positions = map { #warn "  $_ => ".$map_positions->{seqname_to_bacname($_)}."\n";
				    $map_positions->{seqname_to_bacname($_)} || ()
				  } $cluster->get_members
	 ) {
	my $avg_position = sum(@map_positions)/@map_positions;
	{ cluster => $cluster, offset => $avg_position, }
      } else {
	#put unmapped clusters way at the end of the chromosome
	{ cluster => $cluster, offset => 600 }
      }
    } $cluster_set->get_clusters;
}

#takes a versioned bac sequence identifier, changes it into a bac name
sub seqname_to_bacname($) {
  assemble_clone_ident('agi_bac_with_chrom',
		       parse_clone_ident(shift,
					 'versioned_bac_seq'
					)
		      )
}


sub copy_agp_files_to_assembly_dir {
    my ($dir, $chrdata) = @_;

    # copy the AGP files into the complete assembly dir also
    for my $chr_rec (values %$chrdata) {
        my $agp_subdir = File::Spec->catdir( $dir, 'generated_agp' );
        mkpath( $agp_subdir );
        my $agp_target = agp_file( $chr_rec->{chrnum}, 1, $agp_subdir );
        copy( $chr_rec->{agp_file}, $agp_target )
            or die "$! copying $chr_rec->{agp_file} -> $agp_target";
    }
}

#given a mummer output file, build a CXGN::Cluster::ClusterSet out of
#it
sub mummer_to_clusterset {
  my ($mummer_output_file) = @_;

  #here's the cluster set we're building
  my $cluster_set = CXGN::Cluster::ClusterSet->new;
  #  $cluster_set->set_debug(1);

  open (my $f, $mummer_output_file)
    or die "Can't open file $mummer_output_file for reading...";

  my $query = "";
  my $reverse;
  while (<$f>) {
    #print "MUMMER: $_\n";
    chomp;

    if (/^>/) {
      (undef,$query,$reverse) = split;
    } else {
      my ($subject, $query_start, $subject_start, $length) = split;

      $cluster_set->add_match( $query,
			       $subject,
			     );

    }
  }

  return $cluster_set;
}


sub write_assembly_metadata {
    my ($assembly_dir, %metadata) = @_;

    mkpath( $assembly_dir );
    YAML::Any::DumpFile( File::Spec->catfile( $assembly_dir, 'assembly_metadata.yml' ),
                         \%metadata,
                        );
}

sub make_assembly_dir_contig_files {
    my ($assembly_dir,$seqs_index) = @_;

    my $all_contigs_seqio = Bio::SeqIO->new(
        -format => 'fasta',
        -file   => '>'.File::Spec->catfile( $assembly_dir, 'contigs_all.fasta' ),
       );

    my $members_filename = _members_filename( $assembly_dir );
    my %set_seqio;
    my %contigs_seqio;
    foreach my $contig_line ( _slurp_chomp( $members_filename ) ) {
        my ($precluster_name,$contig_name,@members) = split /\t/,$contig_line;
        my $precluster_number = _extract_precluster_number( $precluster_name );
        my ($tag) = $precluster_name =~ /^(\w+)_precluster/;

        # also open a sequence output for each of the sequence sets
        my $set_contigs_seqio =
            $set_seqio{$tag} ||= #< lazily build seqios for each contig set (per chromosome)
                Bio::SeqIO->new(
                    -format => 'fasta',
                    -file   => '>'.File::Spec->catfile( $assembly_dir, "contigs_$tag.fasta" ),
                   );

        if( -d _precluster_dir( $assembly_dir, $tag, $precluster_number ) ) {
            my $precluster_singlets_filename = _precluster_dir( $assembly_dir, $tag, $precluster_number, 'precluster_members.seq.singlets' );
            -s $precluster_singlets_filename
                and die "precluster $precluster_name has nonempty singlets file $precluster_singlets_filename, this should not have anything in it!";

            my $precluster_contigs_filename  = _precluster_dir( $assembly_dir, $tag, $precluster_number, 'precluster_members.seq.contigs' );

            # cache the seqio for each precluster .contigs file, so
            # that as we go through the members file, we will
            # sequentially get the contig seqs out of the contigs
            # file, which should be in the same order as the .members
            # file.  this is a little dangerous
            my $contigs_in = $contigs_seqio{$precluster_contigs_filename}
                ||= Bio::SeqIO->new( -format => 'fasta', -file => $precluster_contigs_filename );

            my $s = $contigs_in->next_seq;
            $s->desc('original_id:'.$s->id);
            $s->id( $contig_name );
            $_->write_seq( $s ) for $all_contigs_seqio, $set_contigs_seqio;

        } else {
            # it's a single sequence in the precluster, just get it from the index
            @members > 1 and die "no phrap .contigs file, but multiple members for precluster $precluster_name, something is wrong";
            my $member = $members[0];
            my $s = $seqs_index->fetch( $member ) or die "what, no $member?";
            $s->desc('original_id:'.$s->id);
            $s->id( $contig_name );
            $_->write_seq( $s ) for $all_contigs_seqio, $set_contigs_seqio;
        }
    }
}

sub make_complete_assembly_tar {
    my ( $dir, $tarfile ) = @_;

    my $bn = basename( $dir );
    local $CWD = File::Spec->catdir( $dir, File::Spec->updir );
    my $cmd = "tar -cf - $bn | gzip -n > $tarfile";
    system $cmd;
    $? and die "$! running '$cmd'\n";
}

sub _gzip_fh {
    # not using IO::Compress::Gzip here because it silently produces
    # corrupt output in some cases (for certain combinations of
    # libraries, I think).  was having this problem with a local-lib
    # running on pipelines.
    my $filename = File::Spec->catfile( @_ );
    open my $fh, "| gzip -nc > '$filename'" or die "$! running gzip";
    return $fh;
}

# slurp a file's lines into an array, chomping each line
sub _slurp_chomp {
    my $file = shift;
    open my $f, $file or die "$! reading $file";
    my @data = <$f>;
    chomp @data;
    return @data;
}

sub _contig_name {
    my ($tag, $ctg_num) = @_;
    return "${tag}_contig$ctg_num";
}
sub _first_number {
    $_[0] =~ /(\d+)/ or die "no number in '$_[0]'";
    return $1;
}
sub _extract_precluster_number {
    $_[0] =~ /precluster\D?(\d+)/ or die "cannot parse '$_[0]'";
    return $1;
}

sub _precluster_name {
    my ($tag,$precluster_number) = @_;
    return "${tag}_precluster${precluster_number}";
}

sub _precluster_dir {
    my ($assembly_dir,$tag,$precluster_number,@additional) = @_;
    return File::Spec->catdir( $assembly_dir, 'precluster_assemblies', _precluster_name( $tag, $precluster_number ), @additional );
}

sub _members_filename {
    my ( $assembly_dir ) = @_;
    return File::Spec->catfile( $assembly_dir, "membership.tab" );
}

# our %seqlens;
# sub sequence_length {
#   my ($sequence_name) = @_;
#   warn "$sequence_name has length $seqlens{$sequence_name}\n";
#   return $seqlens{$sequence_name};
# }
# sub index_sequence_lengths {
#   my ($seqfile) = @_;
#   open my $s,$seqfile or die "$! opening $seqfile";
#   my $i;
#   while(<$s>) {
#     if(s/^>//) {
#       ($i) = split;
#     } else {
#       $seqlens{$i} += length($_) - 1;
#     }
#   }
# }
