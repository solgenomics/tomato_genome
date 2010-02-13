#!/usr/bin/perl
use strict;
use warnings;

#use Data::Dumper;

use CXGN::CDBI::Class::DBI { dbargs => {AutoCommit => 0}}; #this is one xtn
use CXGN::Genomic;

###### CONFIGURATION ######

#order of columns in the file
#change this if your column order changes
my @column_order = qw/genbank_dbgss_id  gss_id  genbank_identifier/;

#pattern that will only be matched by a line with data on it 
#(as opposed to column headings or comments or whatnot)
#change this if your first column is no longer strictly numeric
my $data_pattern = qr'^\s*\d';

##### /CONFIGURATION ######


sub usage {

  die <<EOF;
 Usage:  genomic_gss_genbank_accession_load.pl <file of ID mappings>

 Takes a file of GSS ID to Genbank accessions and loads the Genbank
 accessions into the Genomic database in the GSS table, and updates the MOST
 RECENT Genbank submission record for each GSS it finds.

 Make sure this is what you want (e.g. this is not the right thing if you
 have done more than one  genbank submission for a given gss, but have not
 loaded the response from either of them, this script will not do what you
 need).

 This program right now assumes three columns in the input file:
      dbGSS ID   |   GSS ID   |   GenBank Accession

 That's easy to change though (\@column_order at the top of this script).
EOF

}

$| = 1; #we want to print out the insertion thingies

usage() unless @ARGV;

my $updated = 0;
my $skipped = 0;
print "Inserting Genbank identifiers.\n('.' = 100 identifiers skipped, '|' = 100 identifiers inserted)\n";
while(<>) {
  next unless $_ =~ $data_pattern;

  my %data;
  @data{@column_order} = split /\s+/;
  $data{gss_id} =~ /^\d+$/
    or die "Invalid GSS ID in input file: '$data{gss_id}'";
  $data{genbank_identifier} =~ /^\w+$/
    or die "Invalid Genbank Accession in input file: '$data{genbank_identifier}'";
  $data{genbank_dbgss_id} =~ /^\d+/
    or die "Invalid dbGSS ID in file: '$data{genbank_dbgss_id}'";

  #find the GSS object we're referring to
  my $gss = CXGN::Genomic::GSS->retrieve( $data{gss_id} )
    or die "No GSS found with ID $data{gss_id}\n";

  #find its most recent
  my @subs = $gss->gss_submitted_to_genbank_objects;
  my $recentsub = $subs[-1]
    or die 'No genbank submission record found for GSS with ID '.squote $gss->gss_id;

  if($recentsub->genbank_identifier) {
    sub disc {
      no strict 'refs'; #use symbolic references
      my ($fieldname) = @_;
      sprintf('Discrepancy found: most recent genbank submission record (ID of %s) for GSS with ID %s has fieldname %s, but this input file has %s',
	      squote($recentsub->gss_submitted_to_genbank_id),
	      squote($gss->gss_id),
	      squote($recentsub->$fieldname),
	      squote($data{$fieldname}),
	     );
    }

    #check that the genbank accessions and dbGSS ids match
    $recentsub->genbank_identifier eq $data{genbank_identifier}
      or die disc('genbank_identifier');

    $recentsub->genbank_dbgss_id == $data{genbank_dbgss_id}
      or die disc('genbank_dbgss_id');

    $skipped++;
    print '.' unless $skipped%100;
  } else {
    $recentsub->set(%data);
#    print "updating with\n".Dumper($recentsub);
    $recentsub->update;
    $updated++;
    print '|' unless $updated%100;
  }
}
print "\n$updated submission records updated.\n$skipped were already up to date the database.\n";


CXGN::Genomic::GSS->dbi_commit;


#single-quote a string
sub squote($) {
  return "'".shift()."'";
}

