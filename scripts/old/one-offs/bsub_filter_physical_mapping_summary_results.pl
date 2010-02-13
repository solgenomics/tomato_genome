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

use List::MoreUtils qw/ any /;

#use Data::Dumper;

our %opt;
getopts('',\%opt) or pod2usage(1);

print scalar(<>); #< print header line of input

my @current_rec;
while(<>) {
  chomp;
  my (undef,$clone_name) = 
    my @fields
      = split /\t/, $_;

  no warnings 'uninitialized';
  unless( $clone_name eq $current_rec[0] ) {
    handle( @current_rec ) if @current_rec;
    @current_rec = ($clone_name);
  }
  push @current_rec, \@fields;
}

sub handle {
  my $clone_name = shift;
  my @recs = @_;

  #index the records by method name
  my %idx;
  push @{$idx{$_->[2]}}, $_ for @recs;

  return unless
       $idx{FISH}
    && ( $idx{Overgo} || $idx{BLAST_vs_marker_seqs} )
    && (    any_positions_match( $idx{FISH}, $idx{Overgo} )
	 || any_positions_match( $idx{FISH}, $idx{BLAST_vs_marker_seqs} )
       );

  print join("\t",@$_)."\n" for @recs;
  print "\n";
}

sub any_positions_match {
  my ($recs_1, $recs_2 ) = @_;

  any {
    my $r1_chr = $_->[5];
    any {
      my $r2_chr = $_->[5];
      $r1_chr == $r2_chr
    } @$recs_2
  } @$recs_1;
    
}

__END__

=head1 NAME

bsub_filter_physical_mapping_summary_results.pl - postprocesses the
output of bsub_physical_mapping_summary.pl (see DESCRIPTION), prints
results on STDOUT.

=head1 DESCRIPTION

Looks for BACs that have both FISH data and associations to markers
(via blast or overgo), for which the FISH and marker data agree on
chromosome placement.

=head1 SYNOPSIS

  bsub_filter_physical_mapping_summary_results.pl [options] summary_results_file

  Options:

    none yet

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
