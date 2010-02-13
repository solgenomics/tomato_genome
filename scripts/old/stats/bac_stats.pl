
#!/usr/bin/perl

use strict;
use Bio::SeqIO;
use category;
use Getopt::Std;
use bac;

=head1 bac stats.pl
    
gives statistics about BACs end sequences. 
Addresses the question of annotations differences between bac ends.

options:  -e: evalue cutoff to use for parsing files
         

Parameters:

(1) categories file
(2) Fasta file of BAC sequences
(3) blast results against NR parsed with parseblast.pl
(4) blast results against repeat dataset
   
Output:

Creates the following files:
(1) Name of seq file + .list + evalue containing all the bac names, annotations, classes, coords, etc.
(2) Name of seq file + .frequencies+ evalue , containing the frequencies of the observed classes

=cut

    
##

;

use vars qw( $opt_b $opt_v $opt_e $opt_g $opt_i);

getopts("be:vgi:");

my $categories_file = shift;
my $bac_fasta = shift;
my $intrinsic_file = $opt_i;
my @blast_reports = @ARGV;
my %seq_classes = ();
my %bacs = ();

my $evalue = $opt_e;


if (!$evalue) { 
    $evalue = 1e-10; 
    print STDERR "No evalue specified. Using $evalue.\n";
}
else { 
    print STDERR "Evalue cutoff set to $evalue.\n";
}

if (!@ARGV) { 
    print "Usage: bac_stats.pl classification_file seq_file.fasta blast_report_1 ... blast_report_n\n";
    print "       [note: results of blast_reports n will overwrite annotations from blast_reports n-1]\n";
    exit();
}

print STDERR "Parsing category defs file [$categories_file]...\n";

# build the categories data structure...
my $category_root = category -> new();
$category_root -> set_root();
$category_root -> parse_category_defs_file($categories_file, 1);
$category_root -> set_name("&repeats");


print STDERR "Reading BAC end sequences [$bac_fasta]...\n";

# read bac ids from bac fasta file
# and generate a hash with bac information
#
my $file = Bio::SeqIO -> new ('-format' => 'Fasta', -file => $bac_fasta);
while (my $seqobj = $file -> next_seq() ) {

  # Get line and chop it up in little pieces...
    
    my $accno = $seqobj -> accession_number();
    my $desc  = $seqobj -> desc(); 
    my $seq  = $seqobj -> seq();
    my $read_id   = uc($seqobj -> id());
    
    $seq =~ s/\*//g;
    $seq =~ s/ +//g;
    $seq =~ s/\r//g;
    $seq =~ s/\n//g;
    $seq = uc($seq);
    
    my $gc = $seq =~ tr/GC/GC/;
    $gc = $gc / length($seq);

    my ($bac_id, $direction) = get_bac_name_from_read($read_id, $desc);
    #print "Read: $bac_id, $direction\n";
    if (!exists($bacs{$bac_id})) { 
	$bacs{$bac_id}=bac->new(); 
	$bacs{$bac_id}->set_id($bac_id);

    }
    $bacs{$bac_id}->set_seq($direction, $seq);
    $bacs{$bac_id}->set_has_end($direction); 
    $bacs{$bac_id}->set_GC_content($direction, $gc);
}

$file -> close();

# get the total number of bac end sequences
#
print STDERR "Getting end counts...\n";
my $end_count = 0;
my $T7_count =0;
my $SP6_count =0;
foreach my $k (keys(%bacs)) { 
    $end_count += $bacs{$k}->get_end_seq_count(); 
    if ($bacs{$k}->has_T7()) { $T7_count++; }
    if ($bacs{$k}->has_SP6()) { $SP6_count++; }
}

print STDERR "Total end seq: $end_count. T7: $T7_count. SP6: $SP6_count\n";

# get the intrinsic repeats (from recon analysis)
#
if ($opt_i) { 
    print STDERR "Getting intrinsic repeat information [$intrinsic_file]...\n";
    open (F, "<$intrinsic_file") || die "Can't open $intrinsic_file\n";
    while (<F>) { 
	chomp;
	my ($fam, $id) = split/\t/;
	my ($bac_id, $dir) = get_bac_name_from_read($id);
	$bac_id = uc($bac_id);
	if (!exists($bacs{$bac_id})) { print STDERR "NOT FOUND: BACID: $bac_id\n"; }
	else { $bacs{$bac_id}->set_intrinsic_repeat($dir, 1); }
    }
    close(F);
}

 
foreach my $blast_report (@blast_reports) { 
    print STDERR "Parsing $blast_report...\n";

    parse_blast_file($blast_report, \%bacs, $evalue);

}
   
print STDERR "Parsing annotations...\n";

my %bac_classes = ();
my %bac_end_classes = ();

foreach my $k (keys(%bacs)) { 
    foreach my $end ("T7", "SP6") { 
	my $class = "";
	my $bac_end_class = "";
	#print STDERR "Analyzing ".$bacs{$k}->get_id()." End: $end\n";
	if (exists($bacs{$k})) { 
	    my ($annotation, $evalue, $bitscore) =$bacs{$k}->get_annotation($end);
	    
	    my @parent_terms = ();	    
#	
	    @parent_terms = $category_root->get_parent_terms($annotation);
	    #shift @T7_parents; # remove &ROOT
	    my $class_path = join ";", @parent_terms;
	    
	    #print STDERR "Annotation: $annotation\n";

	    # parent_terms[0] should be the root element...
	    # if there are annotation but no repeat category was detected, it should be a gene
	   
	    if (($annotation ) && (!$parent_terms[1])) { 
		$class = "&gene";
		$bac_end_class = $class;
	    }
	    elsif ($parent_terms[1]) { 
		$class = $class_path;
		if ($parent_terms[1] eq "&TEs") { 
		    $bac_end_class = $parent_terms[1];
		}
		else {
		    $bac_end_class = $parent_terms[1];
		}
	    }
	    elsif ($bacs{$k}->get_intrinsic_repeat($end)) { 
		$class = "&intrinsic";
		$bac_end_class = $class; 
	    }
	    else { 
		$class = "&no_annotation"; 
		$bac_end_class = $class;
	    }
	    
	    # check if there is sequence information for either end. If not, set 
	    # sequence class to &not_sequenced
	    #
	    if (!$bacs{$k}->has_end($end)) { $class="&not_sequenced"; $bac_end_class = $class; }
	    
	    $bacs{$k}->set_end_class($end, $bac_end_class);
	    $bacs{$k}->set_end_class_path($end, $class);
	    
	    $seq_classes{$class_path}++;
	}

	 
    }
    my $T7_class = $bacs{$k}->get_end_class("T7");
    my $SP6_class = $bacs{$k}->get_end_class("SP6");
    
    my $class = join( ", ", sort ($T7_class, $SP6_class));
#	print STDERR "BAC class = $class\n";
    
    $bac_end_classes{$T7_class}++;
    $bac_end_classes{$SP6_class}++;
    
    $bac_classes{$class}++;

    
}
my $total_bac_end_count = 0;

# calculate the theoretical number of classes.
# first, determine the number of bac ends involved
#
foreach my $k (keys %bac_end_classes) { 
    #print "$k\t$bac_end_classes{$k}\n";
    $total_bac_end_count += $bac_end_classes{$k};
}

# then, calculate the theoretical bac classes based on the 
# observed frequencies
#
my %theoretical_bac_classes = ();
foreach my $c1 (sort keys(%bac_end_classes)) { 
    foreach my $c2 (sort keys(%bac_end_classes)) { 
	my ($b1, $b2) =  sort ($c1, $c2);
	
	$theoretical_bac_classes{$b1.", ".$b2} += $bac_end_classes{$b1}*$bac_end_classes{$b2}/($total_bac_end_count)/2;
	#print STDERR "Classes: $b1 [$bac_end_classes{$b1}] $b2 [$bac_end_classes{$b2}] bac end count: $end_count ".$theoretical_bac_classes{$b1.", ".$b2}."\n";
    }
}

# calculate the GC content per class

my %gc_content_per_class = ();

foreach my $k (keys(%bacs)) {
    foreach my $end ("T7", "SP6") { 
	my $bclass = $bacs{$k}->get_end_class($end);
	$gc_content_per_class{$bclass} += $bacs{$k}->get_GC_content($end);
    }
}

# evaluate final results...

my $total =0;
my $total_theoretical = 0;

open(F, ">$bac_fasta.list.$evalue") || die "Can't open $bac_fasta.list.$evalue"; 
 
# print a list of all bacs with their classes.

foreach my $k (keys(%bacs)) { 
    my $has_T7 = $bacs{$k}->has_T7();
    my $has_SP6 = $bacs{$k}->has_SP6();
    my $intrinsic_T7 = "";
    if ($bacs{$k}->get_intrinsic_repeat("T7")) { $intrinsic_T7="***"; }
    my $intrinsic_SP6 = "";
    if ($bacs{$k}->get_intrinsic_repeat("SP6")) { $intrinsic_SP6="***"; }
    print F "$k\t".($bacs{$k}->get_end_class_path("T7"))."\t$intrinsic_T7\t".($bacs{$k}->get_GC_content("T7"))."\t".($bacs{$k}->get_blast_coords("T7"))."\t".($bacs{$k}->get_annotation("T7"))[0]."\t".($bacs{$k}->get_end_class_path("SP6"))."\t$intrinsic_SP6\t".($bacs{$k}->get_GC_content("SP6"))."\t".($bacs{$k}->get_blast_coords("SP6"))."\t".($bacs{$k}->get_annotation ("SP6"))[0]."\n";
    $total++;
}

close (F);

open (F, ">$bac_fasta.frequencies.$evalue") || die "Can't open $bac_fasta.frequencies.$evalue";

print F "Observed frequencies: (class, count, percent, average gc content\n";
foreach my $k (keys %bac_end_classes) { 
    print F "$k\t$bac_end_classes{$k}\t".($bac_end_classes{$k}/$total_bac_end_count)."\t".($gc_content_per_class{$k}/$bac_end_classes{$k})."\n";
    
}

print F "Total: $total_bac_end_count\n\n";

print F "\n\nFrequencies of each class of sequence: \n\n"; 

foreach my $k (sort (keys(%seq_classes))) { 
    print F "$k\t$seq_classes{$k}\t".($seq_classes{$k}/$total_bac_end_count)."\n";
}

print F "\n\n";

$total = 0;
print F "BAC name\tBAC class\tCount\tTheoretical\n";
foreach my $c1 (sort keys(%bac_end_classes)) {
    foreach my $c2 (sort keys(%bac_end_classes)) {
	#my ($class1, $class2) = split /\s*\,\s*/, $k;
	if (! exists($bac_end_classes{$c1})) { print STDERR "$c1 not found...\n"; }
	if (! exists($bac_end_classes{$c2})) { print STDERR "$c2 not found...\n"; }
	
	my ($class1, $class2) =  ($c1, $c2);
	
	my $bac_class = $class1.", ".$class2;
	
	if ($bac_class !~ /not_sequenced/) { 
	    print F "$bac_class\t".$bac_classes{$bac_class}."\t".$theoretical_bac_classes{$bac_class}."\n";
	    $total +=$bac_classes{$bac_class};
	    $total_theoretical += $theoretical_bac_classes{$bac_class};
	}
    }
    
}


# calcualte the GC content for each BAC end and for each BAC etc
# we already calculated GC content for the full BAC end seq
# now lets calculate the GC content for the annotated region on the end
# lets also compare the GC content in either BAC end.
print STDERR "Calculating GC content for annotated regions...\n";
open(GC, ">$bac_fasta.gc_content") || die "Can't open $bac_fasta.gc_content";
foreach my $k (keys(%bacs)) { 
    # extract the annotated sequence
    
    foreach my $dir ("T7", "SP6") { 
	my $subseq = "";
	foreach my $coords ($bacs{$k}->get_blast_coords($dir)) { 
	    my @hits = split /\,/, $coords; 
	    foreach my $h (@hits) { 
		my ($s, $e) = split /\.\./, $coords; 
		#print STDERR "S: $s, End: $e\n";
		
		if ($s > $e) {  #reverse orientation
		    my $s = substr($bacs{$k}->get_seq($dir), $e, ($e+$s-1));
		    my $revcomp = reverse($s);
		    $revcomp =~ tr/ATGC/TACG/;
		    $subseq = $revcomp . $subseq;
		}
		else { 
		    $subseq .= substr($bacs{$k}->get_seq($dir), $s, ($s+$e-1));
		}
	    }
	}
	#print STDERR "Subseq=".(length($subseq))." seq=".(length($bacs{$k}->get_seq($dir)))."\n";
	my $gc_content = $subseq =~ tr/GC/GC/; 
	if (length($subseq)!=0) { 
	    $gc_content = $gc_content/length($subseq); 
	    print GC $bacs{$k}->get_id()."\t$dir\t".($bacs{$k}->get_end_class($dir))."\t$gc_content\n";
       }
	else { print GC $bacs{$k}->get_id()."\t".($bacs{$k}->get_GC_content($dir))."\t(length of seq is 0!)\n"; }
    }
    
}
close(GC);
			     
			
		       




# foreach my $k (keys(%bac_end_classes)) { 
#     $total_theoretical += $bac_end_classes{$k};
# }

# correct total theoretical classes, &not_sequenced, &not_sequenced BACs were
# not in the original dataset
# 
$total_theoretical = $total_theoretical - $theoretical_bac_classes{"&not_sequenced, &not_sequenced"};

print F "Total BACs analyzed: $total. Theoretical total: $total_theoretical\n";

close(F);

print STDERR "Done. Output is in $bac_fasta.frequencies.$evalue and $bac_fasta.list.$evalue\n";

sub parse_blast_file {
    # parses a blast file that has been parsed using parseblast.pl 
    # (contains id in col 0, evalue in col 10 an annotation in col 12).
    my $file = shift;
    my $bacs_ref = shift;
    my $EVALUE_CUTOFF = shift;

    if (!$EVALUE_CUTOFF) { 
	$EVALUE_CUTOFF = 1e-10; 
	print STDERR "No evalue cutoff was specified. Using default [$EVALUE_CUTOFF]\n";
    }
    

    open (NR, "<$file") || die "Can't open $file\n";
    
    my $previous_id = "";
    my $best_match_id = "";
    my $previous_direction = "";
    my $previous_bac_id = "";
    my @blast_coords = ();
    while (<NR>) { 
	chomp;
	my @f = split /\t/;
	my ($id, $match_id, $q_start, $q_end, $evalue, $bitscore, $annotation) = ($f[0], $f[1], $f[6], $f[7], $f[10], $f[11], $f[13]);
	$id = uc($id);
	my ($bac_id, $direction) = get_bac_name_from_read($id);       
	if ($previous_id ne $id) {
	    
	    # first line for new id, contains the most significant hit
	    
	    if ($previous_bac_id) { 
		#$$bacs_ref{$previous_bac_id}->clear_blast_coords($previous_direction);
		foreach my $b (@blast_coords) { 
		    my ($s, $e) = split /\s*\,\s*/, $b;
		    $$bacs_ref{$previous_bac_id}->add_blast_coords($previous_direction, $s, $e);
		}
	    }
	    @blast_coords =();
	    if ($evalue < $EVALUE_CUTOFF) { 
		#print "BAC_ID: $bac_id, $annotation, $evalue\n";
#
		
		my($old_annotation, $old_evalue, $old_bitscore) = $$bacs_ref{$bac_id}->get_annotation($direction);

		if ($old_bitscore < $bitscore) { 
		    
		    $$bacs_ref{$bac_id}->set_annotation($direction, "$match_id [$evalue, $bitscore] $annotation", $evalue, $bitscore);
		    $best_match_id = $match_id;
		}
	    }
	    else { 
		$best_match_id = "";
	    }
	}
	if ($match_id eq $best_match_id) {
	    push @blast_coords, "$q_start,$q_end";
	}
	$previous_id=$id;
	$previous_direction = $direction;
	$previous_bac_id=$bac_id;
    }
}

sub get_bac_name_from_read { 
    my $read_name = shift;
    my $bac_id;
    my $direction;
     if ($read_name =~ /(.*)(SP6|T7)/i) { 
	$bac_id = $1;
	$direction = $2;
	#print STDERR "SPLIT $read_name into $bac_id and $direction...\n";
	return ($bac_id, $direction);
    }
    if ($read_name =~ /(.*)(f|r)$/i) { 
	$bac_id = $1;
	$direction = $2;
	if ($direction =~ /f/i) { $direction ="T7"; }
	if ($direction =~ /r/i) { $direction ="SP6"; }
	return ($bac_id, $direction);
    }
    
    else {
	#print STDERR "Couldn't split $read_name...\n";
	return $read_name; 
    }
}
