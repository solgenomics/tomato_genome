#!/usr/bin/perl

=head1 NAME

generate_agp_file.pl - a quick script to generate AGP files from the sequencing data provided to SGN

=head1 DESCRIPTION

generate_agp_file [options] bac_fasta_file

options:

=over 4

=item -F

Don't perform the splitting of the sequence file into different chromosome - use data in the agp_tmp directory instead (from a previous run of this script).

=item -C

Don't run mummer again. Use the result files in the agp_tmp directory instead (from a previous run of this script).

=back

The script will generate an agp file for each chromosome, and name them sgn_chrXX.agp. They should be copied to the /data/prod/public/tomato_genome/agp/ directory, from where the mapviewer will try to read them if no country-provided file is available. Note that the files generated are not quite legal agp files, but sufficient for the mapviewer to show the data.

Cheers!

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut

use strict;
use File::Spec;
use Getopt::Std;
use Bio::SeqIO;
use CXGN::DB::Connection;
use CXGN::Cview::MapFactory;
use CXGN::Cluster;

our ($opt_F, $opt_C);
getopts('FC');

my $bac_fasta = shift;
my $mummer = "/usr/bin/mummer";
my $min_overlap = 1000;


my $dbhost = "scopolamine";
my $dbname = "cxgn_tmp";

# separate out the BACs according to chromosome
#
my $TEMPDIR = "agp_tmp";
mkdir ($TEMPDIR);

my @chr = map { sprintf "%02d", $_; } (1..12);


my %out_file = ();
my %out = ();
if ($opt_F) { 
    print STDERR "-F option: not generating the files.\n";
}
else { 
    my $in = Bio::SeqIO->new( -format=>"fasta",
			      -file => $bac_fasta
			      );

    print STDERR "Partitioning file into chromosomes...\n";


    foreach my $chrf (@chr) { 
	print "Generating output file chr $chrf...\r";
	$out_file{$chrf} = File::Spec->catfile("$TEMPDIR", "BACS-$chrf.fasta");
	$out{$chrf} = Bio::SeqIO->new( -format => "fasta",
				       -file => ">".$out_file{$chrf},
				       );
	
    }
    
    while (my $bac = $in->next_seq()) { 
	
	print STDERR "Processing ".($bac->id())."\r";
	if ($bac->id() =~ /C(\d{2}).*/i) { 
	    my $chr = $1;
	    #print STDERR "Chromosome : $chr\n";
	    $out{$chr}->write_seq($bac);
	}
    }


    # close files
    foreach my $k (keys(%out)) { 
	$out{$k}->close();
    }
    $in->close();
}
# run mummer on each chromosome file with itself...
#
if ($opt_C) { 
    print STDERR "Not running mummer... attempting to use old files.\n"; 
}
else { 
    print "Running mummer...\n\n";
    foreach my $chrf (@chr) {
	print "Running mummer on chr $chrf ($out_file{$chrf})...\n";
	system ("$mummer -mum -b -l $min_overlap -n -F $out_file{$chrf} $out_file{$chrf} > $TEMPDIR/mummer_out_$chrf");
    }
}
# construct the clusters
# ...
print STDERR "Constructing the clusters...\n";
my %cluster_set = ();

local($!);
$!=1;
foreach my $chrf (@chr) { 
    $cluster_set{$chrf} = CXGN::Cluster::ClusterSet->new();
    $cluster_set{$chrf}->set_debug(1);
    my $F;
    print "Analyzing chr $chrf...\n";
    open ($F, "<$TEMPDIR/mummer_out_$chrf") || die "Can't open file mummer_out_$chrf for reading...";
    my $query = "";
    while (<$F>) { 
	chomp;
	
	print ".";
	my $reverse = 0;
	my $subject = "";
	my $query_start = 0;
	my $subject_start = 0;
	my $length = 0;
	if (/^>/) { 
	    if (/REVERSE/i) { 
		$reverse = 1;
	    }
	    if (/^> (.*)/) { 
		$query = $1; 
	    } 
	    if ($query =~ /(.*)\s+Reverse/) { 
		$query = $1;
	    }
	    #print STDERR "QUERY: $query\n"; 
	}
	else { 
	    (undef, $subject, $query_start, $subject_start, $length) = split /\s+/;
	    #print STDERR "READING: $subject, $query_start, $subject_start, $length\n";
	    #print STDERR "QUERY: $query. SUBJECT: $subject\n";

	    if ($subject =~ /(.*)\.\d+$/) { 
		$subject = $1;
	    }
	    if ($query =~ /(.*)\.\d+$/) { 
		$query = $1;
	    }

	    $cluster_set{$chrf}->add_match($query, $subject, $query_start, ($query_start+$length), $subject_start, $subject_start+$length);

	}	    

	
    }
    foreach my $cluster ($cluster_set{$chrf}->get_clusters()) { 
	print "\nCLUSTER ".($cluster->get_unique_key())."\n";
	$cluster->calculate_contig_coords();
	$cluster->get_contig_coords();
	foreach my $member ($cluster->get_members()) { 
	    print "$member\n";
	}
	
    }

}


#exit(-1);
# get marker associations from SGN
print STDERR "Pulling in marker information...\n";
my $dbh = CXGN::DB::Connection->new( { "dbhost"=> $dbhost,
                                       "dbname" => $dbname,
				       "dbuser" => "web_usr",
				       "dbpass" => "sol\@ley!" }
				     );
			
#make CXGN::Genomic::Clone happy	    
$ENV{DBHOST}=$dbhost;
$ENV{DBNAME}=$dbname;
$ENV{DBSCHEMA}="sgn";
$ENV{DBBRANCH}="devel";

my $map_factory = CXGN::Cview::MapFactory->new($dbh);
my $map= $map_factory->create( {map_version_id => "p9"});
my %bacs = ();
foreach my $chrf (@chr) { 
    print STDERR "Get map information for chromosome $chrf (".(int($chrf)).")\n";
    my $chr = $map->get_chromosome(int($chrf));
    print STDERR "Chr $chrf has ".scalar($chr->get_markers)." markers on it.\n";
    foreach my $bac ($chr->get_markers()) { 
	
	if ($bac->isa("CXGN::Cview::Marker::Physical")) { 
	    $bacs{$bac->get_marker_name()} = $bac->get_offset();
	}
	#print "BAC: ".($bac->get_marker_name())."\n";
	
    }
    
    my %mapped_clusters = ();
    my %cluster_positions = ();

    print "Getting clusters for $chrf...\n";
    foreach my $cluster ($cluster_set{$chrf}->get_clusters()) { 
	foreach my $member ($cluster->get_members()) { 
	    
	    if (exists($bacs{$member})) { 
		print "BAC $member is associated to chromosome $chrf and has position $bacs{$member}.\n";
		$mapped_clusters{$cluster->get_unique_key()}++;
		$cluster_positions{$bacs{$member}}= $cluster;
	    }
	    else { 
		print "No association found for $member.\n";
	    }
	}
    }

    my $mapped_count = 0;
    my $total_count = 0;
    foreach my $cluster ($cluster_set{$chrf}->get_clusters()) { 
	$total_count++;
	if (exists($mapped_clusters{$cluster->get_unique_key()})) { 
	    #print STDERR $cluster->get_unique_key()." is mapped\n"; 
	    $mapped_count++;
	    
	}
	else { 
	    #print STDERR $cluster->get_unique_key()." is NOT mapped\n";
	}

    
    }

    print "Total mapped: $mapped_count out of $total_count\n";

    # we assume that 1cM is on average about 750kB of sequence (Wing, Zhang, and Tanksley, 1994).
    #
    my $correspondence = 200000;

    # output an approximation of the AGP file...
    #    
    my $out;
    my $out_path = File::Spec->catfile($TEMPDIR, "sgn_chr$chrf.agp");
    open ($out, ">$out_path") || die "Can't open $out_path for writing";
    
    print STDERR "Writing agp information to $out_path...\n";
    
    my $count = 1;
    foreach my $offset (sort { $a <=> $b } (keys(%cluster_positions))) { 
	my $base = $offset * $correspondence;
	if ($mapped_clusters{$cluster_positions{$offset}->get_unique_key()} > 1) { 
	    # determine order... TBD
	}

	foreach my $info ($cluster_positions{$offset}->get_contig_coords()) { 
	    my ($name, $contig_start, $contig_end) = @$info;
	    my $global_start = $base + $contig_start -1;
	    my $global_end   = $base + $contig_end -1;
	    print $out "S.lycopersicum-chr$chrf\t$global_start\t$global_end\t$count\tF\t$name\n";
	}
	$count++;
    }

    close ( $out);
}



