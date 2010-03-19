#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Std;
use POSIX; #for numerical stuff
use UNIVERSAL qw/isa/;

use DBI;
#DBI->trace(2,'/tmp/dbitrace');

use CXGN::DB::Connection;
CXGN::DB::Connection->verbose(0);

use CXGN::Page::FormattingHelpers qw/html_break_string/;

use CXGN::Genomic::Library;
use CXGN::Genomic::GenbankSubmission;
use CXGN::Genomic::Search::GSS;

################ CONFIGURATION VARS ################

#default for min hqi len
my $default_min_hqi_len = 200;
#number of GSS files to put in a single output directory.  GSS files
#are output in directories named for the ranges their GSS IDs fall
#into.
my $gss_files_per_directory = 10_000;

####################################################
sub usage {
  die <<EOU;

genomic_gss_genbank_submit.pl - for each library specified,
  check which sequences still need to be genbank submitted,
  create genbank submission package for them, and make a
  gss_genbank_submission record for them.

Usage:

  genomic_gss_genbank_submit.pl [ options ] -L <libname>,... -d <dirname>
                                            -t 'publication title'
                                            -c 'contact name'

  note: -L, -d, -t, and -c options are all REQUIRED

Options:

  -d  <dirname>
      Directory in which to put the submission files.

  -L  <shortname>,<shortname>,...
      List of library shortnames to process for possible genbank submissions.
      If the special name 'ALL' is given, makes submissions for all libraries
      in the database.

  -t  <title>
      Title of the publication these submissions will appear in.

  -c  <contact name>
      Name of the contact person to give to Genbank.

  -u  <username>
      Record this as being done by a different user name.  Default is the
      real user running the script.

  -x  Dry run.  Only generate files, do not touch the database,
      do not record this as being sent to GenBank.

  -l  <# base pairs>
      Minimum high-quality region length for a sequence to be considered
      worth submission.  Default $default_min_hqi_len.

  -f  Do not enforce the 'flags == 0' requirement on GSS seqs that are output,
      that is, output GSS objects regardless of the contents of their gss.flags
      field.

EOU
}

sub progprint(@); #prints unbuffered progress messages

my %opt;
getopts('d:L:c:t:l:xf', \%opt) || usage();

#check that mandatory options are set
$_ or usage() foreach (@opt{qw/ d L c t /});

#default min hqi len
$opt{l} = $default_min_hqi_len unless defined $opt{l};

#check if the target directory exists.  if not, create it
my $targetdir = $opt{d};
-w $targetdir
  || mkdir($targetdir)
  || die "Could not create/write to target directory $targetdir\n";

#make sure we properly start our transaction
#NOTE: this script is one big transaction
CXGN::CDBI::Class::DBI->db_Main->dbh_param( AutoCommit => 1 );
CXGN::CDBI::Class::DBI->db_Main->begin_work;

#get the list of library objects from which we'll be submitting
my @libs = do {
  if( $opt{L} =~ /^all$/i ) {
    CXGN::Genomic::Library->retrieve_all;
  } else {
    my @shortnames = split ',',$opt{L};
    @shortnames or usage(); #must have a list here
    map { CXGN::Genomic::Library->search( shortname=>$_ ) } @shortnames;
  }
};

#make a row in the database to record this submission batch
my $submission_object =
  CXGN::Genomic::GenbankSubmission-> create({ date_generated        => strftime("%Y-%m-%d %H:%M:%S",localtime),
					      submitted_by          => scalar(getpwent),
					    });

foreach my $lib (@libs) {
  progprint "Making Genbank submission for library ".$lib->shortname.".\n";
  my $libid = $lib->library_id;
  my $libdir = File::Spec->catdir($targetdir,$lib->shortname);
  mkdir $libdir; #may already be there, that's fine

  #make a Publication, Library, and Contact file for this submission
  make_static_genbank_files($opt{c}, $opt{t}, $lib, $libdir);

  #do a search for the GSS objects we will be dumping
  my $matching_gss = do {
    progprint "\tQuerying GSS objects to dump...";
    my $gss_search = CXGN::Genomic::Search::GSS->new;
    #should probably be able to buffer
    #200_000k GSS objects at a time
    $gss_search->page_size(200_000);
    #$gss_search->debug(1);
    my $gssq = $gss_search->new_query;
    $gssq->library_id("=$libid");
    #  $gssq->gss_id("< 20");
    $gssq->trimmed_length(" >= $opt{l}") if $opt{l} && $opt{l} > 0;
    $gssq->flags('=0') unless $opt{f};
    $gssq->status_not_set('vec_unk | contam_unk');
    $gssq->needs_genbank_submit;
#    $gssq->gss_id('<200000');
    my $gssr = $gss_search->do_search($gssq);
    $gssr->autopage($gssq,$gss_search);
    progprint "done\n";
    $gssr;
  };

  #dump each GSS we need to dump
  my $gsscount = 0;
  my $gss_filename = File::Spec->catfile( $libdir, $lib->shortname.'.gss' );
  open(my $gss_fh, ">$gss_filename")
    or die "Could not open '$gss_filename' for writing: $!\n";
  progprint "\tDumping ".$matching_gss->total_results." GSS objects into $gss_filename...\nDumped 0...";

  my $gss;
  while($gss = $matching_gss->next_result) {
    $gsscount++;
    $gsscount % 200 or progprint "$gsscount...";

    write_gss_genbank_record($gss, $lib, $opt{c}, $opt{t}, $gss_fh );

    #record this as submitted
    $submission_object->add_to_gss_submitted_to_genbank({ gss_id => $gss });
  }
  progprint "\n\tDumped total of $gsscount GSS objects from lib ".$lib->shortname.".\n";
}

#will have died and rolled back before here if not committed
unless( $opt{x} ) {
  CXGN::CDBI::Class::DBI->db_Main->commit;
}

=head1 FUNCTIONS

=head2 make_static_genbank_files

  args: contact name string,
        publication title string,
        CXGN::Genomic::Library object,
        destination directory
  ret : 1 on success, 0 on failure

  submission file formats are specified at:
  http://www.ncbi.nlm.nih.gov/dbGSS/how_to_submit.html

=cut

sub make_static_genbank_files {
  my ($contactname,$pubtitle,$lib,$targetdir) = @_;
  progprint "\tMaking static contact, library, and publication files in $targetdir.\n";

  my $retval = 1;
  $retval &&= make_contact_genbank_file($contactname,"$targetdir/Cont");
  $retval &&= make_pub_genbank_file($pubtitle,"$targetdir/Pub");
  $retval &&= make_lib_genbank_file($lib,"$targetdir/Lib");
  return $retval;
}

=head2 make_contact_genbank_file

  Desc: make a genbank-format Contact file, with the passed
        contact name filled in (fill in the rest yourself after the
        script runs)
  Args: contact name, filename for new template contact file
  Ret : 0 on failure, 1 on success

=cut

sub make_contact_genbank_file {
  my ($contactname,$targetfile) = @_;

  #make a Contact file
  # TYPE:   Entry type - must be "Cont" for contact entries.
  #         **Obligatory field**.
  # NAME:   Name of person providing the GSS sequence
  #         **Obligatory field**.
  # FAX:    Fax number as string of digits.
  # TEL:    Telephone number as string of digits.
  # EMAIL:  E-mail address
  # LAB:    Laboratory
  # INST:   Institution name
  # ADDR:   Address string
  # ||
# warn "making contact file $targetfile\n";
  open(CON,">$targetfile")
    or return 0;
  print CON genbank_format_file( [ type  => 'Cont'                ],
				 [ name  => $contactname          ],
				 [ fax   => 'insert fax number w/area code'],
				 [ tel   => 'insert telephone number w/area code'],
				 [ email => 'rmb32@cornell.edu'   ],
				 [ lab   => 'Dr. Steven Talksley' ],
				 [ inst  => 'Cornell University'  ],
				 [ addr  => '252A Emerson Hall, Cornell University, Ithaca, NY 14850 USA' ],
			       );
  close CON;
  return 1;
}

=head2 make_pub_genbank_file

  Desc:
  Args: filename for new pub file
  Ret : 0 on failure, 1 on success

=cut

sub make_pub_genbank_file {
  my ($pubtitle,$targetfile) =@_;
  #make a Publication file
  # TYPE:    Entry type - must be "Pub" for publication entries.
  #          **Obligatory field**.
  # MEDUID:  Medline unique identifier.
  # 	 Not obligatory, include if you know it.
  # TITLE:   Title of article.
  #          **Obligatory field**.
  # 	 Begin on line below tag, use multiple lines if needed
  # AUTHORS: Author name, format:  Name,I.I.; Name2,I.I.; Name3,I.I.
  #          **Obligatory field**.
  # 	 Begin on line below tag, use multiple lines if needed
  # JOURNAL: Journal name
  # VOLUME:  Volume number
  # SUPPL:   Supplement number
  # ISSUE:   Issue number
  # I_SUPPL: Issue supplement number
  # PAGES:   Page, format:   123-9
  # YEAR:    Year of publication.
  #          **Obligatory field**.
  # STATUS:  Status field.1=unpublished, 2=submitted, 3=in press,
  # 	 4=published
  #          **Obligatory field**.
  # ||

  open(PUB,">$targetfile")
    or return 0;
  print PUB genbank_format_file( [ type     => 'Pub'                                  ],
				 [ meduid   => ''                                     ],
				 [ title    => $pubtitle ,                'multiline' ],
				 [ authors  => 'Mueller,L.; Tanksley,S.', 'multiline' ],
				 [ journal  => 'Fake Journal Name'                    ],
				 [ volume   => 'Fake volume number'                   ],
				 [ suppl    => 'supplement number'                    ],
				 [ issue    => 'issue'                                ],
				 [ i_suppl  => 'issue supplement number'              ],
				 [ pages    => '41-9'                                 ],
				 [ year     => '2005'                                 ],
				 [ status   => 1                                      ],
			       );
  close PUB;
  return 1;
}

=head2 make_lib_genbank_file

  Desc: Create a Library file for a submission to NCBI Genbank's dbGSS
  Args: Library object, filename for the library file to make
  Ret : 0 on failure, 1 on success

=cut

sub make_lib_genbank_file {
  #make a Library file
  # TYPE:      Entry type - must be "Lib" for library entries.
  #            **Obligatory field**.
  # NAME:      Name of library.
  #            **Obligatory field**.
  # ORGANISM:  Organism from which library prepared.
  # STRAIN:    Organism strain
  # CULTIVAR:  Plant cultivar
  # ISOLATE:   Individual isolate from which the sequence was obtained
  # SEX:       Sex of organism (female, male, hermaphrodite)
  # ORGAN:     Organ name
  # TISSUE:    Tissue type
  # CELL_TYPE: Cell type
  # CELL_LINE: Name of cell line
  # STAGE:     Developmental stage
  # HOST:      Laboratory host
  # VECTOR:    Name of vector
  # V_TYPE:    Type of vector (Cosmid, Phage, Plasmid, YAC, other)
  # RE_1:      Restriction enzyme at site1 of vector
  # RE_2:      Restriction enzyme at site2 of vector
  # DESCR:     Description of library preparation methods,
  # 	   vector, etc.
  #            This field starts on the line below the DESCR: tag.
  # ||

  my ($lib,$targetfile) = @_;
  open(LIB,">$targetfile")
    or return 0;
  my ($accession,$organism,$common_name) = $lib->accession_name;
  $organism = join(' ', ( split /\s+/, $organism )[0..1] );
  print LIB genbank_format_file( [ type     => 'Lib'                             ],
				 [ name     => $lib->name                        ],
				 [ organism => $organism                         ],
				 [ cultivar => $common_name                      ],
				 [ host     => $lib->cloning_host                ],
				 [ vector   => $lib->vector                      ],
				 [ v_type   => uc($lib->clone_type_object->name) ],
				 [ re_1     => $lib->rs1                         ],
				 $lib->rs2 ?
				([ re_2     => $lib->rs2                         ])
                                              : (),
			       );
  close LIB;
  return 1;
}


=head2 write_gss_genbank_record

  Desc: dump the given L<CXGN::Genomic::GSS> object into a GSS file
        suitable for submitting to genbank
  Args: GSS object, Library object, contact name string,
        publication title string,
        filehandle to append to
  Ret : 0 on failure, 1 on success

=cut

sub write_gss_genbank_record {
  my ($gss,$lib,$contactname,$pubtitle,$target) = @_;
  #make GSS sequence files for each sequence
  # TYPE:          Entry type - must be "GSS" for GSS entries.
  #                **Obligatory field**
  # STATUS:        Status of GSS entry - "New" or "Update".
  #                **Obligatory field**
  # CONT_NAME:     Name of contact
  #                Must be identical string to the contact entry
  # 	       **Obligatory field**
  # CITATION:      Journal citation
  # 	       Must be identical string to the publication title
  #                Begins on line below tag.
  #                Use continuation lines if needed.
  # 	       **Obligatory field**
  # LIBRARY:       Library name
  #                Must be identical string to library name entry.
  # 	       **Obligatory field**
  # GSS#:          GSS name or number assigned by contact lab. For GSS entry
  #                updates, this is the string we match on.
  # 	       **Obligatory field**
  # GDB#:          Genome Database accession number
  # GDB_DSEG:      Genome Database Dsegment number
  # CLONE:         Clone number/name
  # SOURCE:        Source providing clone, e.g., ATCC
  # SOURCE_DNA:    Source identity number for the clone as pure DNA
  # SOURCE_INHOST: Source identity number for the clone stored in the host
  # OTHER_GSS:     Other GSSs on this clone.
  # DBNAME:        Database name for cross-reference to another
  # 	       database
  # DBXREF:        Database cross-reference accession
  # PCR_F:         Forward PCR primer sequence
  # PCR_B:         Backward PCR primer sequence
  # INSERT:        Insert length (in bases)
  # ERROR:         Estimated error in insert length (bases)
  # PLATE:         Plate number or code
  # ROW:           Row number or letter
  # COLUMN:        Column number or letter
  # SEQ_PRIMER:    Sequencing primer description or sequence
  # P_END:         Which end sequenced, e.g., 5'
  # HIQUAL_START:  Base position of start of high-quality sequence
  #                (default = 1)
  # HIQUAL_STOP:   Base position of last base of high-quality
  # 	       sequence
  # DNA_TYPE:      Genomic (default), cDNA, Viral, Synthetic, Other
  # CLASS:         Class of sequencing method, e.g., BAC ends,
  # 	       YAC ends, exon-trapped
  # 	       **Obligatory field**
  # PUBLIC:        Date of public release
  # 	       Leave blank for immediate release.
  # 	       **Obligatory field**
  #                Format:   MM/DD/YYYY
  # PUT_ID:        Putative identification of sequence by submitter
  # COMMENT:       Comments about GSS.
  #                Text starts on line below COMMENT: tag.
  # SEQUENCE:      Sequence string.
  #                Text starts on line below SEQUENCE: tag.
  # 	       **Obligatory field**
  # ||	

  my $chromat = $gss->chromat_object;
  my $clone = $chromat->clone_object;
  my $qc = $gss->qc_report_object;

  my @other_gss;
  foreach my $chr ( $clone->chromat_objects ) {
    my $g = $chr->latest_gss_object;
    push @other_gss, $g unless !$g || $g->gss_id == $gss->gss_id;
  }
  @other_gss = map {$_->gss_id} @other_gss;

  my ($hqs,$hqlen,$hqe);
  if($qc) {
    $hqs = $qc->hqi_start;
    $hqlen = $qc->hqi_length;
    $hqe = $hqs+$hqlen-1;
  }

  my $seq = html_break_string($gss->seq || '',60,"\n");

  print $target genbank_format_file( [ type          => 'GSS'                    ],
				     [ status        => 'New'                    ],
				     [ cont_name     => $contactname             ],
				     [ citation      => $pubtitle   , 'multiline'],
				     [ library       => $lib->name               ],
				     [ 'gss#'        => $gss->gss_id             ],
				     [ clone         => $clone->clone_name       ],
#				     [ source        => 'fake source name'       ],
# 				     [ source_dna    => ''                       ],
# 				     [ source_inhost => ''                       ],
				     [ other_gss     => join(',', @other_gss)    ],
				     [ insert        => $clone->estimated_length ],
#				     [ error         => ''                       ],
				     [ plate         => $clone->platenum         ],
				     [ row           => $clone->wellrow          ],
				     [ column        => $clone->wellcol          ],
				     [ seq_primer    => $chromat->primer         ],
				     [ p_end         => $chromat->direction      ],
				     #   dbGSS flappily starts indexing at 1
				     [ hiqual_start  => $hqs + 1                 ],
				     [ hiqual_stop   => $hqe + 1                 ],
				     [ dna_type      => 'Genomic'                ],
				     [ class         => $chromat->read_class_object->class_name ],
#				     [ public        => ''                       ],
#				     [ comment       => ''                       ],
				     [ sequence      => $seq        , 'multiline'],
				   );
  return 1;
}

=head2 genbank_format_records

  Desc:
  Args: record field name, string of data to be put in the field,
        flag telling whether multiline or not (perl-evaluated to true/false)
  Ret : a string containing the properly formatted genbank submission
        field, suitable for printing to a file or stdout or whatever

=cut

sub genbank_format_record {
  my ($fieldname,$data,$multiline) = @_;
  my $retstring;
  $data ||= '';


  $retstring = uc($fieldname).': ';
  $retstring .= $multiline ? "\n$data" : $data;
  $retstring .= "\n" unless $retstring =~ /\n$/;
  return $retstring;
}

=head2 genbank_format_file

  Desc: Iterate genbank_format_record over arguments, close the file with a
        '||' on a line by itself, to make a valid genbank-submission-format
        file.
  Args: Any number of array refs, whose contents will be fed as arguments
        directly to genbank_format_record
  Ret:  string containing full text of a valid genbank submission file

=cut

sub genbank_format_file {
  my @fieldrecs = @_;
  my $retstr;

  foreach my $fieldrec (@fieldrecs) {
    $retstr .= genbank_format_record(@$fieldrec);
  }
  $retstr .= "||\n";
  $retstr;
}


sub progprint(@) {
  local $| = 1;
  print @_;
}
