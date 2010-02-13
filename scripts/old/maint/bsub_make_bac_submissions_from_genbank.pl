#!/usr/bin/env perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;

use English;
use Carp;
use FindBin;

use Getopt::Std;
use Pod::Usage;

use Path::Class;

use File::Temp;
use File::Spec::Functions;

use Bio::DB::GenBank;
use Bio::SeqIO;

use CXGN::Debug;
use CXGN::Genomic::Clone;
use CXGN::Genomic::CloneIdentifiers qw/parse_clone_ident/;

use Data::Dumper;

our %opt;
getopts('',\%opt) or pod2usage(1);
# my $d = CXGN::Debug->new;

# if( $opt{x} ) {
#     $d->set_debug(1);
# }

my $gb_recs  = Bio::DB::GenBank->new(-retrievaltype => 'tempfile');

while( <> ) {
    chomp;
    next unless /\S/; #< skip any blank lines

    my @f = split;
    my $input_clone_name = shift @f;
    my $gbacc = shift @f;

    eval {

        my $clone = CXGN::Genomic::Clone->retrieve_from_clone_name( $input_clone_name )
            or die "cannot find clone for '$input_clone_name'";
        my $clone_name = $clone->clone_name;
        $clone->chromosome_num
            or die "clone '$input_clone_name / $clone_name' has no chromosome number in the BAC registry";


        # fetch the full genbank seq
        my $gb_seq = $gb_recs->get_Seq_by_acc( $gbacc )
            or die "could not fetch accession '$gbacc'";

        # make sure the input clone name matches the genbank record's clone name
        my $clone_2 = _infer_clone_obj_from_gb_richseq( $gb_seq )
            or die "could not infer clone from gb record:\n".Dumper($gb_seq);
        $clone->clone_id == $clone_2->clone_id
            or die "input claims $gbacc is clone '$clone_name', but Genbank record has clone '".$clone_2->clone_name."'";

        # make a BAC submission tarball for it
        my $clone_name_with_chrom = $clone->clone_name_with_chromosome;

        my $tempdir = File::Temp->newdir;
        my $bac_dir = catdir( $tempdir, $clone_name_with_chrom );

        mkdir $bac_dir or die "could not make dir '$bac_dir'";
        # write the genbank accession
        file($bac_dir,'gbacc.txt')->openw->print($gb_seq->accession."\n");

        $gb_seq->primary_id( $clone_name_with_chrom ); #< change the primary seq ID to the clone name with chr
        # write the sequence file
        Bio::SeqIO->new( -format => 'fasta',
                         -file => '>'.catfile($bac_dir,"$clone_name_with_chrom.seq")
                       )
                  ->write_seq( $gb_seq );
        my $tarball_name = "$clone_name_with_chrom.tar.gz"; #< will be made in the current dir
        system 'tar', -C => $tempdir, -czf => $tarball_name, "$clone_name_with_chrom/"
            and die "failed: tar -C $tempdir -czf $tarball_name $clone_name_with_chrom/ ($!/$?)";
    };
    if( $EVAL_ERROR ) {
        warn "failed ( $input_clone_name  $gbacc ):\n$EVAL_ERROR";
    }
}

sub _infer_clone_obj_from_gb_richseq {
    my $seq = shift;
    my $clone_name = _infer_clone_name_from_gb_richseq($seq)
        or return;

    require CXGN::Genomic::Clone;
    return CXGN::Genomic::Clone->retrieve_from_clone_name( $clone_name );
}

sub _infer_clone_name_from_gb_richseq {
    my $seq = shift;

    my @tests =
        (
         sub { shift->desc =~ /clone\s+([\w\-]+)/ && $1 },
         sub { my $a = shift->annotation(); my ($c) = $a->get_Annotations('clone'); $c && $c->value },
        );

    foreach my $t (@tests) {
        my $name = $t->($seq);
        return $name if $name && CXGN::Genomic::CloneIdentifiers::parse_clone_ident($name);
    }

    return;
}

__END__

=head1 NAME

bsub_make_bac_submissions_from_genbank.pl <query_list>

=head1 SYNOPSIS

  bsub_import_clone_seqs_from_genbank.pl [options] gb_queries_file

  NOTE: can also read genbank queries from stdin

  Options:

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
