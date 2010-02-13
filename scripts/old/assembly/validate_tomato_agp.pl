#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use CXGN::Tools::List qw/ min max /;
use File::Temp qw/tempfile/;
use File::Basename;

use Bio::PrimarySeq;

use CXGN::BioTools::AGP qw/agp_parse/;
use CXGN::DB::Connection;
CXGN::DB::Connection->verbose(0); #< shut up
use CXGN::Tools::Run;

######## DEFAULTS ############

our $default_identity_threshold = 99;
our $default_min_overlap_length = 100;

##############################


sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script  agp_file   agp_file   ...

  Given a list of tomato genome AGP files, check each of them for
  syntax (like the more general validate_agp.pl script), but also do a
  homology check on all overlapping sequence regions.

  Options:

   -i <pct>
     minimum percent identity for overlapping regions
     default: $default_identity_threshold

  -m <length>
    minimum overlap length for checking to be performed
    default: $default_min_overlap_length

EOU
}
sub HELP_MESSAGE {usage()}


our %opt;
getopts('i',\%opt) or usage();

my $identity_threshold = defined $opt{i} ? $opt{i} : $default_identity_threshold;
my $min_overlap_length = defined $opt{m} ? $opt{m} : $default_min_overlap_length;

foreach my $file (@ARGV) {

  my $bn = basename($file);

  my @errors;
  my $error_flag;
  my $error = sub { warn "$bn:".shift()."\n";  $error_flag = 1 };

  # parse the AGP file, validate its syntax and identifiers
  my $agp_lines = agp_parse($file, validate_syntax => 1, validate_identifiers => 1 );
  unless( $agp_lines ) {
    warn "$file failed syntax/identifier validation.\n";
    next;
  }

  # now go through the lines two at a time, and do bl2seq on each
  # overlap, and warning if any regions have percent identity less
  # than the threshold
  my $previous_line;
  foreach my $current_line (@$agp_lines) {
    next if $current_line->{comment}; #< skip comments

    if( $current_line->{ident} ) {
      # fetch the sequence for this line
      my $seq = $current_line->{seq} ||= fetch_seq( $current_line->{ident} );

      # validate the lengths of this component
      my $seqlen = $seq->length;
      $current_line->{cend} <= $seqlen
	or $error->("$current_line->{linenum}: component end $current_line->{cend} is beyond end of sequence ($seqlen bases)");

      # see if it overlaps with the previous component.  if so,
      # validate the overlapping sequence
      if( $previous_line && $previous_line->{ident} ) {
	# figure out the length of the sequence overlap
	my $overlap_length = do {
	  my $prev_o_ind = $previous_line->{ostart} - $previous_line->{cstart} + 1;
	  my $curr_o_ind =  $current_line->{ostart} -  $current_line->{cstart} + 1;
	  my $curr_ol = $previous_line->{oend} - $curr_o_ind + 1;
	  my $prev_ol = $prev_o_ind + $previous_line->{seq}->length - 1 - $current_line->{ostart} + 1;
	  #warn "line: $current_line->{linenum} prev: $prev_ol curr: $curr_ol\n";
	  $curr_ol + $prev_ol
	};

	$overlap_length = min( $overlap_length, $current_line->{seq}->length, $previous_line->{seq}->length );

	# now compare the overlapping sequences
	if( $overlap_length >= $min_overlap_length ) {
	  # extract each sequence
	  my $prev_ov_seq = $previous_line->{seq}->trunc( $previous_line->{seq}->length - $overlap_length + 1,
							  $previous_line->{seq}->length );
	  $prev_ov_seq = $prev_ov_seq->revcom if $previous_line->{orient} eq '-';
	  my $curr_ov_seq =  $current_line->{seq}->trunc( 1, $overlap_length );

	  my $id_pct = find_id_pct( $prev_ov_seq, $curr_ov_seq );
	  #warn "$current_line->{linenum}: $previous_line->{ident} to $current_line->{ident} : $overlap_length bases, $id_pct% identical\n";
	  unless( $id_pct >= $identity_threshold ) {
	    my $id1 = $previous_line->{ident};
	    my $id2 = $current_line->{ident};
	    $error->("$current_line->{linenum}: overlapping sequences ($id1 and $id2, $overlap_length bases of overlap) are only $id_pct% identical");
	  }

	}
      }
    }

    $previous_line = $current_line;
  }

}


# given a sequence identifier, return a Bio::PrimarySeqI containing
# its sequence, or die if fetch failed


sub fetch_seq {
  my ($seqname) = @_;

  # look up the sequence
  our $dbh ||= CXGN::DB::Connection->new;
  my $seqs = $dbh->selectall_arrayref(<<EOQ,undef,$seqname);
select name,residues from public.feature where name = ?
EOQ

  # check that we got exactly 1 sequence back
  @$seqs > 1
    and die "multiple rows in public.feature table found with name '$seqname'\n";
  @$seqs
    or die "no sequences found with name '$seqname'\n";

  return Bio::PrimarySeq->new( -id => $seqs->[0][0], -seq => $seqs->[0][1] );
}


# given two Bio::PrimarySeqI, find their percent identity
sub find_id_pct {
  my ($seq1, $seq2) = @_;

  my ($t1_fh,$t1_f) = tempfile( File::Spec->catfile( File::Spec->tmpdir, 'validate-tomato-agp-seq1-XXXXXX'), UNLINK => 0 );
  my ($t2_fh,$t2_f) = tempfile( File::Spec->catfile( File::Spec->tmpdir, 'validate-tomato-agp-seq1-XXXXXX'), UNLINK => 0 );

  print $t1_fh '>'.$seq1->id."\n".$seq1->seq."\n";
  print $t2_fh '>'.$seq2->id."\n".$seq2->seq."\n";

  #warn '>'.$seq1->id."\n".$seq1->seq."\n";
  #warn '>'.$seq2->id."\n".$seq2->seq."\n";

  close $t1_fh;
  close $t2_fh;

  #warn "running bl2seq...\n";
  my $bl2 = CXGN::Tools::Run->run( 'bl2seq',
				   -i => $t1_f,
				   -j => $t2_f,
				   -D => 1,
				   -p => 'blastn',
				 );

  unlink $t1_f, $t2_f;

  open my $o, $bl2->out_file or die "$! reading ".$bl2->out_file;
  my $max_pct_id = 0;
  while( my $l = <$o> ) {
    next if $l =~ /^\s*#/;
    my (undef,undef,$pct) = split /\s+/,$l;
    $max_pct_id = $pct if $pct > $max_pct_id;
  }
  close $o;

  return $max_pct_id;
}
