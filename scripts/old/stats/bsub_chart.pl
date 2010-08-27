#!/usr/bin/perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;

use English;
use Carp;
use FindBin;

use Fatal qw/ open /;
use Getopt::Std;
use Pod::Usage;
use File::Copy;

use DateTime;
use Path::Class;

use Chart::Lines;

use CXGN::TomatoGenome::BACPublish qw/ parse_filename /;
use CXGN::Tools::List qw/ flatten /;
use CXGN::TomatoGenome::Config;

use Data::Dumper;

my %opt;
getopts('r:cl:m:v',\%opt) or pod2usage(1);

my $listing_source = do {
    if( $opt{l} ) {
        $opt{l}
    } else {
        my $remote_ssh = $opt{r} ? "ssh $opt{r}" : '';
        my $pubdir = dir( map CXGN::TomatoGenome::Config->new->get_conf($_), 'ftpsite_root', 'bac_publish_subdir' );
        unless( $remote_ssh || $opt{l} ) {
            -d $pubdir or die "Cannot find bac publishing dir '$pubdir'.  Please check ftpsite_root and bac_publish_subdir config variables";
            -r $pubdir or die "Cannot read from dir '$pubdir'";
        }

        "$remote_ssh ls -lt --time-style=+\%s $pubdir/chr*/*/C*.tar.gz |";
    }
};

# list the BAC tarballs and parse their dates, counting them into bins
#warn "listing bacs...\n";
my %counts;
open( my $ls, $listing_source );
my $min_date;
my $max_date;
while(my $line = <$ls>) {
    chomp $line;
    my @f = split /\s+/, $line, 7;
    unless( $f[6] && parse_filename( $f[6] ) ) {
        warn "invalid line: $line\n";
        next;
    }
    my $date = DateTime->from_epoch( epoch => $f[5] )
        or die "cannnot convert from epoch $f[5]";
    $min_date = $date unless $min_date && $min_date < $date;
    $max_date = $date unless $max_date && $max_date > $date;
    my $bin = datetime_to_bin_name($date);
    $counts{$bin}++;
}
close $ls;

$min_date && $max_date or die "no data, cannot generate graph";

#warn "calculating bins...\n";
my @all_bins = all_bins( $min_date, $max_date );
my @data = map $counts{$_->epoch} || 0, @all_bins;

if( $opt{c} ) {
    my $c = 0;
    @data = map $c += $_, @data;
}

if( $opt{v} ) {
    for( my $i = 0; $i < @all_bins; $i++ ) {
        print STDERR $all_bins[$i]->mdy.'  '.$data[$i]."\n";
    }
}

if( $opt{m} ) {
    my ($m,$d,$y) = split m!/!,$opt{m},3;
    my $mindate = DateTime->new( year => $y, month => $m, day => $d )
        or die "cannot parse -m '$opt{m}'";
    while( @all_bins && $all_bins[0] < $mindate ) {
        shift @all_bins;
        shift @data;
    }
    @all_bins or die "the given -m '$opt{m}' excludes all data, cannot graph.\n";
}

my @bin_labels = map { $_->strftime('%b/%y') } @all_bins;

my $chart = Chart::Lines->new( 800,600 );
$chart->set(
            ### configuration settings for this chart
            title => 'Uploaded BACs by Month',
            skip_x_ticks => 31,
            x_grid_lines => 1,
            precision => 0,
            legend => 'none',
           );

$chart->png( shift(@ARGV) || \*STDOUT, [\@bin_labels, \@data] );


########## SUBROUTINES ##############
sub all_bins {
    my ($min_date,$max_date) = @_;
    $min_date = $min_date->clone;
    my $one_day = DateTime::Duration->new( days => 1 );
    my %all_bins;
    while( $min_date <= $max_date ) {
        my $bin = datetime_to_bin_name( $min_date );
        $all_bins{ $bin } = 1;
        $min_date->add_duration( $one_day );
    }

    return
        map DateTime->from_epoch( epoch => $_ ),
        sort {$a <=> $b} keys %all_bins;
}

# minimum bin size is one day
sub datetime_to_bin_name {
    my $d = shift;

    # this is, for example, 2009-03
    my $date = $d->clone;
    #if( $date->day > 15 ) {
    #    $date->subtract_duration( DateTime::Duration->new( months => 1) );
    #}
    $date->set(
               #day => 1,
               hour => 0,
               minute => 0,
               second => 0,
              );

    #warn "$d -> $date\n";
    return $date->epoch;
}

__END__

=head1 NAME

bsub_chart.pl - graph BAC uploads

=head1 DESCRIPTION

Uses the file modification times of BAC submission tarballs in the
main ftp site dir to graph BAC submissions over time.  Produces an
800x600 PNG graphic file, either to the given file name or on stdout.

=head1 SYNOPSIS

  bsub_chart.pl [options] target_filename.png
  bsub_chart.pl [options] > target_filename.png

  Options:

    -r hostname
       read BAC directories on the given host over ssh

    -c
      if passed, make bac numbers cumulative

    -l file
      if passed, read the bac file listing from the given file.
      conflicts with r

    -m date
       MM/DD/YYYY string that sets a lower bound for data to graph

    -v if passed, print the data points being graphed to stderr

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
