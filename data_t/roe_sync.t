use strict;
use warnings;

use Test::More;

use LWP::Simple;

use CXGN::Genomic::CloneIdentifiers qw/ parse_clone_ident assemble_clone_ident /;


use_ok('CXGN::TomatoGenome::CmdLine');
use_ok('CXGN::TomatoGenome::CmdLine::Command::roe_sync');

my $roe_html = get('http://www.genome.ou.edu/tomato_table.html');
open my $roe_f, '<',  \$roe_html or die "cannot open string ref as filehandle";
while( my $line = <$roe_f> ) {
    next unless $line =~ /BAC/;
    my ($bac) = $line =~ /\s([a-z]{3}-\S+)\s/;
    my ($acc) = $line =~ /(AC\d+)/;

    my $p = parse_clone_ident( $bac );
    ok( $p, "parsed roe clone ident '$bac'" );
    my $agi_name = assemble_clone_ident( agi_bac => $p );
    ok( $agi_name, 'assembled agi name' );

    #check that it is on the SGN FTP
    like( ftp_data( $agi_name )->{accession}, qr/^$acc\.\d+$/, 'Roe and FTP accessions match' );

    #fetch the page on SGN and check that its sequence and genbank accession are linked
    my $detail_url = "http://solgenomics.net/search/quick_search.pl?term=$bac";
    my $detail_page = get($detail_url);
    like( $detail_page, qr!>\[Download fasta\]</a>!, "detail page for $bac has sequence available ($detail_url)" );
}

done_testing;
exit;

######## SUBROUTINES ########

my %ftp_data;
sub ftp_data {
    my ($agi_name) = @_;
    _get_published_bacs() unless %ftp_data;
    return $ftp_data{$agi_name};
}

sub _get_published_bacs {
    my $accs = get('ftp://ftp.sgn.cornell.edu/tomato_genome/bacs/curr/bacs_accessions.txt');
    open my $accs_f, '<', \$accs or die;
    while( <$accs_f> ) {
        chomp;
        my ($bac, $acc) = split;
        my $p = parse_clone_ident( $bac );
        ok( $p, "parsed SGN FTP clone ident '$bac'" );
        my $agi_name = assemble_clone_ident( agi_bac => $p );
        ok( $agi_name, "assembled agi name for SGN FTP clone ident '$bac'" );
        $ftp_data{$agi_name} = {
            'with_chrom' => $bac,
            'accession' => $acc,
        };
    }
}
