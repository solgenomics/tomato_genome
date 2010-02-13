#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;
use DateTime;
use DateTime::Format::Strptime;

use CXGN::DB::Connection;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script

  This script does something.

  Options:

   -g [ day | week | month ]
      time granularity for binning data

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('g:',\%opt) or usage();
$opt{g} ||= 'day';
$opt{g} = lc $opt{g};

my $dbh = CXGN::DB::Connection->new;

my $increment_duration = DateTime::Duration->new( $opt{g}.'s' => 1);

my %ip;
my %ip_and_complete;

print qq|"Date","In Progress","Complete"\n|;

my %bac_stats;
my %stat_bacs;

my $status_log = $dbh->selectall_arrayref("select bac_id,status,timestamp from sgn_people.bac_status_log order by timestamp");

my $ip = 0;
my $ip_and_complete = 0;
my $dt_fmt = DateTime::Format::Strptime->new( pattern => '%F %T', time_zone => 'America/New_York', locale => 'en_US' );

#decorate the returned rows with their parsed timestamps
foreach my $row (@$status_log) {
  $row->[3] = $dt_fmt->parse_datetime($row->[2])
    or die "cannot parse timestamp '$row->[2]'";
}

#this sub takes a row and updates the bac counts with it
my $process_row = sub {
  my ($row) = @_;
  my ($bac_id,$status,$time) = @$row;
  $status ||= 'none';
  my $dt = $dt_fmt->parse_datetime($time)
    or die "could not parse timestamp '$time'";

  $stat_bacs{ $bac_stats{$bac_id} }-- if $bac_stats{$bac_id};
  $bac_stats{$bac_id} = $status;
  $stat_bacs{$status} ||= 0;
  $stat_bacs{$status}++;
};

my $begin_ts = $status_log->[0]->[2];
$begin_ts =~ s/\d+:\d+:\d+\.\d+$/0:0:0/;

my $curr_time = $dt_fmt->parse_datetime($begin_ts);

#now march through time in set increments, processing rows as we go, and printing out the total at each time increment
while( @$status_log ) {
  while( @$status_log &&  DateTime->compare($status_log->[0]->[3],$curr_time) < 1 ) {
    $process_row->(shift @$status_log);
  }

  print join ',',$curr_time->ymd('-'),map $_||0,$stat_bacs{'in_progress'},$stat_bacs{'complete'};
  print "\n";

  $curr_time->add( $increment_duration );
}
