
use strict;
use Bio::SeqIO;

my $file = shift;

open (F, "<$file") || die "Can't open $file";
my @file = <F>;
close(F);

my $count=0;
foreach my $line (@file) {
    chomp($line);
    my ($bac, $genbank_accession) = split /\t/, $line;
    my $chr;
    my $bacname;
    
    print STDERR "Processing $bac, $genbank_accession\n";

    if ($bac =~ /^C(\d+)([A-Z]+)(\d+[A-Z]\d+)$/i) { 
	$chr=$1;
    }
    else { 
	print "Not a valid BAC name! try again!\n";
	exit(-1);
    }
    if ($bac =~ /^(C\d+[A-Z]+)(\d{3}[A-Z]\d+)$/i) { 
	$bacname = $1."0".$2;
    }
    else { 
	$bacname=$bac;
    }
    my $update = "";
    if ($genbank_accession) { 
	$update = " -a $genbank_accession "; 
    }

    my $fasta = "/data/shared/tomato_genome/bacs/chr$chr/$bac/$bac.seq";
    if (-e "$fasta.screen") { 
	#print STDERR "Taking the screened sequence.\n";
	#$fasta = $fasta.".screen";
    }

    my $out = "/data/shared/tomato_genome/bacs/genbank_submission/Cornell:$bacname.sqn";

    print STDERR "Input file: $fasta\n";
    print STDERR "Output file: $out\n";

    $count++;
    system("fa2htgs -i $fasta -t ~/sequin/cornell_template.sqn -u T -p 3 -g Cornell -b 0 -s $bacname -c $bacname -h $chr -C HBa $update -n \"Lycopersicon esculentum\" -o $out");

}


print STDERR "Ran $count times.\n";
