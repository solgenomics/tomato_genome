#!/usr/bin/env perl
use strict;
use warnings;
#use autodie ':all';

use English;
use Carp;
use FindBin;
use File::Path;
use List::MoreUtils qw/ any /;

use Path::Class;

use Memoize;

use Getopt::Std;
use Pod::Usage;

use Bio::DB::GenBank;
use Bio::SeqIO;

use CXGN::Genomic::Clone;
use CXGN::Genomic::CloneIdentifiers qw/ parse_clone_ident assemble_clone_ident /;
use CXGN::Debug;
our $d = CXGN::Debug->new;

use CXGN::PotatoGenome::Config;
use CXGN::PotatoGenome::FileRepository;

use CXGN::DB::DBICFactory;

use Data::Dumper;

our %opt;
getopts('R:vU',\%opt) or pod2usage(1);
sub vsay(@) {
    print @_,"\n" if $opt{v} || $d->get_debug;
}

my $cfg = CXGN::PotatoGenome::Config->load_locked;

$opt{R} ||= $cfg->{repository_path};

# connect to our DB
use Bio::Chado::Schema;
Bio::Chado::Schema->load_classes('MyCloneFeature');

# log in with the name of the user running this unless another one has been specified
$ENV{DBUSER} = '' unless defined $ENV{DBUSER};
my $chado = CXGN::DB::DBICFactory->open_schema('Bio::Chado::Schema', config => $cfg, search_path => ['genomic','public'], );
#$chado->storage->dbh_do(sub { my ($sp) = $_[1]->selectrow_array('SHOW search_path'); warn "search path is $sp\n"} );

load_potato_bacs_into_chado( $chado );
update_potato_bacs_ftp_repos( $chado, $opt{R} );

exit;

############ SUBROUTINES ###########

sub update_potato_bacs_ftp_repos {
    my ( $chado, $repos_dir ) = @_;
    $repos_dir = Path::Class::Dir->new( $repos_dir ) unless ref $repos_dir;

    # open/make our potato BAC repos
    my $repos = CXGN::PotatoGenome::FileRepository
                  ->new( basedir => $repos_dir,
                         create => 1,
                       );

    -w $repos_dir->stringify
        or die "dir '$repos_dir' is not writable";

    # go through all our individual clone seqs and make seq tempfiles
    # for those that are in the DB but not yet in the file repos, and
    # return publish operations for them
    { my @publish_ops = assemble_publish_operations_for_single_clone_seqs( $chado, $repos );
      $repos->publish( @publish_ops );
    }


    ### now generate aggregate files for all clone seqs, and Sol ID -> Genbank mapping
    { my @publish_ops;

      # all clone seqs
      my $aggregate_vf = $repos->get_vf( class => 'AllCloneSequences',
                                         format => 'fasta.gz',
                                       )
          or die 'could not find AllCloneSequences vf';

      # genbank <-> Sol ID mapping
      my $mapping_vf = $repos->get_vf( class => 'AccessionMapping',
                                       format => 'txt',
                                     )
          or die 'could not find AccessionMapping vf';

      my $new_mapping   = File::Temp->new;
      my $new_aggregate = File::Temp->new;
      $new_aggregate->close;
      open $new_aggregate, "| gzip -nc > $new_aggregate" or die $!;
      foreach my $seqfile ( sort $repos->search_files( class => 'SingleCloneSequence' ) ) {
          my $slurp = $seqfile->slurp;
          my ($seqname,$genbank_accession) = $slurp =~ /^\s*>\s*(\S+).+genbank_accession:(\S+)/
              or die "could not extract seqname and genbank acc from $seqfile";
          $new_aggregate->print( $seqfile->slurp );
          $new_mapping->print("$seqname\t$genbank_accession\n");
      }
      $_->close for $new_mapping, $new_aggregate;

      push @publish_ops,
          $aggregate_vf->publish_new_version( $new_aggregate ),
          $mapping_vf->publish_new_version($new_mapping);


      # now publish a new version of the aggregate files
      $repos->publish( @publish_ops );
    }
}

#extract the seq identifier from the file (a Path::Class::File obj)
sub _extract_single_defline {
    my $file = shift;
    my $fh = $file->openr(); #< openr dies on failure
    while (<$fh>) {
        return $_ if /^\s*>/;
    }
    die "cannot extract defline from '$file'";
}

sub assemble_publish_operations_for_single_clone_seqs {
    my ( $chado, $repos ) = @_;

    vsay "searching for existing files in repos...";

    # find which seqs already have seq files
    my %existing_seqs_in_repos =
      # and for each of those, find the seq_id inside the file (should
      # probably be on the first line)
        map {
            my $seq_vf = $_;
            my $seq_file = $_->current_file;

            #extract the seq identifier from the file
            my ($seq_id) = _extract_single_defline($seq_file) =~ /^\s*>\s*(\S+)/;

            # return $seq_id => 1 to make the hash of seq ids that are already present
            $seq_id
                or die "'$seq_file' does not seem to be a valid fasta file!";

            #check that the file in the repos looks plausibly good
            if( $seq_file->stat->size > 5_000 ) {
                $seq_id => $seq_vf
            } else {
                #ignore it and try to dump it again if it does not look good
                warn "file too small, trying to dump again: $seq_file\n";
                ()
            }

        }
          # filtered for vfs that have current files
        grep $_->current_file, 
          # find all single-clone-sequence vfs
        $repos->search_vfs( class => 'SingleCloneSequence' );

    # this is the list of publishing operations we'll be assembling
    my @publish_ops;

    # get Potato org and BAC_clone cvterm
    my $potato = potato_organism($chado);
    #my $bac_cvterm = bac_cvterm($chado);

    vsay "searching for clone seqs in DB...";

    # go through all the current clone sequences and
    my $clone_seqs = $chado->resultset('MyCloneFeature')
                           ->search_related('feature',
                                            { organism_id => $potato->organism_id }
                                           );

    while( my $clone_seq = $clone_seqs->next ) {
        my $seq_id = $clone_seq->name;
        # skip if this seq has already been published
        if( !$opt{U} && $existing_seqs_in_repos{$seq_id} ) {
            vsay "$seq_id: already in repos, skipping";
            next;
        }
        vsay "$seq_id: adding to repository";

        my $clone = CXGN::Genomic::Clone->retrieve_from_clone_name($seq_id)
            or die "failed to retrieve clone for '$seq_id'";

        # otherwise, make publishing ops for any old versions that are
        # published for this clone
        my @rm_ops =
            map $_->publish_remove,
            grep $_->current_file,
            $repos->search_vfs( class => 'SingleCloneSequence',
                                clone => $clone,
                              );
        push @publish_ops, @rm_ops;

        my $project_country = $chado->resultset('Cv::Cvterm')
                                    ->search({ name => 'project_country'})
                                    ->search_related( 'featureprops',
                                                      { feature_id => $clone_seq->feature_id },
                                                      { order_by => 'type_id DESC' },
                                                    )
                                    ->first;
        unless( $project_country ) {
            die "no project_country featureprop found for clone seq '"
                .$clone_seq->name
                ."', feature_id "
                .$clone_seq->feature_id;
        }

        # find the vf where this seq should go
        my $vf = $repos->get_vf( class => 'SingleCloneSequence',
                                 sequence_name => $clone_seq->name,
                                 project => $project_country->value,
                                 format => 'fasta',
                               );

        # make a new seq file for this one
        my $new_seq = File::Temp->new;
        Bio::SeqIO->new( -fh => $new_seq, -format => 'fasta')
                  ->write_seq( $clone_seq );
        $new_seq->close;

        push @publish_ops, $vf->publish_new_version( $new_seq );
    }

    return @publish_ops;
}

# this sub queries genbank for all potato BACs it has and creates or
# updates them in our chado DB as necessary so that our chado reflects
# what is in genbank

sub load_potato_bacs_into_chado {
    my $chado = shift;

    # set up two genbank handles, one that gets records with no seqs, and
    # one that gets just fasta seqs
    my $gb_recs  = Bio::DB::GenBank->new(-retrievaltype => 'tempfile' ,
                                         #-format => 'Fasta',
                                         #don't download any of the sequences
                                         -seq_start => 1, -seq_stop  => 1,
                                        );
    my $gb_fasta = Bio::DB::GenBank->new( -retrievaltype => 'tempfile',
                                          -format => 'Fasta',
                                        );

    vsay "querying GenBank for potato bacs...";

    my $potato_bac_recs = $gb_recs->get_Stream_by_query( Bio::DB::Query::GenBank->new
                                                         ( -db => 'nucleotide',
                                                           #-query  => 'AC236750',
                                                           -query => 'RHPOTKEY or POTGEN or RH89-039-16',
                                                         )
                                                       );
    my $count = 0;
    while ( my $seq = $potato_bac_recs->next_seq ) {

        # find the record for the BAC associated with this seq
        my $clone = find_bac_from_gb_richseq( $seq );
        unless( $clone ) {
            warn "cannot find clone record for ".$seq->accession_number.", skipping.\n";
            next;
        }

        my $project = infer_project_from_gb_richseq( $seq )
	    or next;
        my $chromosome = infer_chromosome_from_gb_richseq( $seq );
        my $phase   = infer_htgs_phase_from_gb_richseq( $seq );

        my $upstream_accession = $seq->accession_number;
        my $upstream_version = $seq->version;

        # find most recent accession for this BAC
        if ( my $current_accession = $clone->genbank_accession( $chado ) ) {

            # compare to this one.  if up to date, next.
            my ($our_version) = $current_accession =~ /\.(\d+)$/
                or die "error parsing genbank accession $current_accession";

            if ( $upstream_version > $our_version ) {
                my $current_seq_name = $clone->latest_sequence_name;
                vsay $clone->clone_name.": current stored version $current_seq_name / $current_accession, loading new genbank seq version $upstream_accession.$upstream_version";
            } else {
                vsay $clone->clone_name.": current stored version $current_accession is >= upstream version $upstream_accession.$upstream_version, skipping";
                next;
            }
        } else {
            vsay $clone->clone_name.": no current sequence, loading upstream $upstream_accession.$upstream_version";
        }

        # fetch the full genbank seq for this one
        my $gb_seq = $gb_fasta->get_Seq_by_gi( $seq->primary_id )
            or die 'could not fetch GI '.$seq->primary_id.":\n".Dumper $seq;

        unless( $gb_seq->length > 5_000 && $gb_seq->seq !~ /[<>]/ ) {
            warn "failed to fetch seq for ".$seq->display_id." (gi ".$seq->primary_id."), skipping.\n";
            next;
        }

        # figure out the new sol-style sequence version for this bac
        my $parsed_ident = parse_clone_ident($clone->clone_name)
            or die "could not parse clone identifier ".$clone->clone_name;
        my $new_sol_seq_name = assemble_clone_ident( versioned_bac_seq_no_chrom =>
                                                     {
                                                      %$parsed_ident,
                                                      version => ($clone->latest_sequence_version || 0) + 1,
                                                     }
                                                   );

        # load it into the DB with the proper accession
        # and update the public.clone_feature table
        $chado->txn_do( sub {

                            my $potato = potato_organism($chado);

                            my $bac_cvterm = bac_cvterm($chado);


                            my $gbacc = $seq->accession_number.'.'.$seq->version;
                            my $gi = $seq->primary_id;

                            # make a feature for it in the feature table
                            my $new_feature =
                                $potato->create_related('features',
                                                        { name => $new_sol_seq_name,
                                                          uniquename => $new_sol_seq_name,
                                                          type_id  => $bac_cvterm->cvterm_id,
                                                          residues => $gb_seq->seq,
                                                          seqlen => $gb_seq->length,
                                                        },
                                                       );
                            $new_feature->create_featureprops({ htgs_phase => $phase,
                                                                finished_seq => ($phase == 3 ? 1 : 0),
                                                                project_country => $project,
                                                                description => "genbank_gi:$gi genbank_accession:$gbacc sequenced_by:$project htgs_phase:$phase chromosome:$chromosome",
                                                                chromosome => $chromosome,
                                                              },
                                                              { autocreate => 1 },
                                                             );


                            # add a dbxref for its genbank accession
                            my $gbacc_dbx = $chado->resultset('General::Db')
                                                  ->find_or_create({ name => 'DB:GenBank_Accession'})
                                                  ->find_or_create_related('dbxrefs',
                                                                           { accession => $gbacc,
                                                                             version   => $seq->version,
                                                                           }
                                                                          );
                            $new_feature->add_to_secondary_dbxrefs( $gbacc_dbx );

                            # add a dbxref for its genbank GI
                            my $gi_dbx = $chado->resultset('General::Db')
                                               ->find_or_create({ name => 'DB:GenBank_GI'})
                                               ->find_or_create_related('dbxrefs',
                                                                        { accession => $gi,
                                                                         }
                                                                       );
                            $new_feature->add_to_secondary_dbxrefs( $gi_dbx );


                            # manually update the clone_feature table
                            $chado->storage->dbh_do(sub {
                                                        my ($s,$dbh) = @_;
                                                        $dbh->do('delete from clone_feature where clone_id = ?',
                                                                 undef,
                                                                 $clone->clone_id,
                                                                );
                                                        $dbh->do('insert into clone_feature (feature_id,clone_id) values (?,?)',
                                                                 undef,
                                                                 $new_feature->feature_id,
                                                                 $clone->clone_id,
                                                                );
                                                    });

                            # don't actually do anything to the db if debug is on
                            $chado->storage->txn_rollback if $d->get_debug;
                        });


        #print Dumper $seq;
        #die unless $clone_name;
        #     print join "\t", $seq->id, $clone_name || 'UNKNOWN', $seq->accession_number, $seq->version;
        #     print "\n";
        $count++;

        # for each seq, we need to update:
        #  - mapping between gb accession and bac name
        #  - bac sequence

    }
    vsay "loaded $count seqs\n";

}


# get the potato organism from chado
sub potato_organism {
    shift()
        ->resultset('Organism::Organism')
        ->find_or_create({ genus => 'Solanum', species => 'Solanum tuberosum'});
}

# get the cvterm for BAC clones from chado
sub bac_cvterm {
    shift()
        ->resultset('Cv::Cv')
        ->search({'me.name' => 'sequence'})
        ->search_related(cvterms => {'cvterms.name' => 'BAC_clone'})
        ->first;
}


sub find_bac_from_gb_richseq {
    my $seq = shift;
    my $clone_name = find_clone_name($seq);

    return CXGN::Genomic::Clone->retrieve_from_clone_name( $clone_name );
}

sub find_clone_name {
    my $seq = shift;
    my @tests =
        (
         sub { shift->desc =~ /clone\s+(\w+)/ && $1 },
         sub { my $a = shift->annotation(); my ($c) = $a->get_Annotations('clone'); $c && $c->value },
        );

    foreach my $t (@tests) {
        my $name = $t->($seq);
        return $name if $name && parse_clone_ident($name);
    }
    return;

}


# # use wget to update the local mirror of bacregistry.potatogenome.net,
# # making sure it preserves modtimes
# #update_local_bacregistry_mirror( $opt{R}_);
# sub update_local_bacregistry_mirror {

#     # for each fasta file in the mirror,
#       # get its file modtime
#       # compare to our local file's modtime

#     # get a list of BACs that are new or changed

#     my $dest_dir = shift;

#     -d $dest_dir
#         or mkpath( $opt{R} )
#         or die "could not create dir '$opt{R}'";

#     my @subdirs = qw( CN IE IN NL NZ PE PL RU US wgs );
#     for my $subdir ( @subdirs ) {
#         system wget =>
#             -P => $dest_dir,# put files at the path given to this script with -R
#             -A => '.fasta', # only download .fasta files
#             -q  =>          # be quiet
#             #-nv =>         # do not print complete verbose progress diags
#             -nH =>          # no hostname prefix in dest dir
#             #-c =>
#             -r =>           # recursive download
#             -N =>           # timestamping
#             -w => 3,        # wait seconds between requests
#             '--random-wait' => # wait between 0.5 and 1.5 times as many seconds as -w above
#             -U => 'SGN Potato Genome Mirror Tool (sgn-feedback@sgn.cornell.edu)',
#             -l => 1,
#             '--no-parent' => #< do not ascend to parent dirs
#             #"http://sgn.localhost.localdomain/data/test_potmirror/",
#             "http://bacregistry.potatogenome.net/pgscreg/downpub/$subdir/fasta";
#     }
# }

sub infer_project_from_gb_richseq {
    my $seq = shift;

    my @project_signatures =
        (
         [ IE => qr/Ireland/i ],
         [ CN => qr/China/, qr/Qu,D/, qr/Du,Y/, qr/He,J/, qr/Zhang,Z/,qr/Huang,S/],
         [ IN => qr/India/],
         [ NL => qr/Netherlands/],
         [ NZ => qr/New Zealand/],
         [ PE => qr/Peru/],
         [ PL => qr/Poland/],
         [ RU => qr/Russia/],
         [ US => qr/United States of America|\bUSA\b/],
         [ AR => qr/Argentina/],
        );

    #match one of the project signatures from the above array
    foreach my $p (@project_signatures) {
        my ( $project_name, @signatures ) = @$p;
        foreach my $sig (@signatures) {
            if( ref $sig eq 'CODE' ) {
                return $project_name if $sig->($seq);
            }
            elsif( ref $sig eq 'Regexp' ) {
                my @refs = $seq->annotation->get_Annotations('reference');
                my $text = join "\n", map {
		    my $ref = $_;
		    map "$_: ".($ref->$_ || '<none>'),
		        qw| title authors location consortium |
		  } @refs;
                return $project_name if $text =~ $sig;
            }
            else {
                die "no handler for '".ref($sig)."' signatures";
            }
        }

    }

    require Data::Dumper;
    warn "WARNING: could not infer project name for sequence:\n".Data::Dumper::Dumper $seq;
    return;
}


sub infer_htgs_phase_from_gb_richseq {
    my $seq = shift;

    my @keywords = map $_->value, $seq->annotation->get_Annotations('keyword','keywords');
    return 1 if any { $_ eq 'HTGS_PHASE1' } @keywords;
    return 2 if any { $_ eq 'HTGS_PHASE2' } @keywords;

    # if no HTGS_PHASE1 or 2, assume phase 3
    return 3;
}

sub infer_chromosome_from_gb_richseq {
    my $seq = shift;

    if( $seq->desc =~ /chromosome[\s\W]+(\d+)/ ) {
        return $1;
    }

    return 'unknown';
}

# on-the-spot DBIC CloneFeature object.  this needs to be somewhere
# other than just here, but for now here it is
BEGIN {
    package Bio::Chado::Schema::MyCloneFeature;
    use base 'DBIx::Class';
    __PACKAGE__->load_components("Core");
    __PACKAGE__->table("clone_feature");
    __PACKAGE__->add_columns(
                             "clone_feature_id",
                             {
                              data_type => "integer",
                              default_value => "nextval('clone_feature_clone_feature_id_seq'::regclass)",
                              is_auto_increment => 1,
                              is_nullable => 0,
                              size => 4,
                             },
                             "feature_id",
                             {
                              data_type => "integer",
                              default_value => undef,
                              is_foreign_key => 1,
                              is_nullable => 1,
                              size => 4,
                             },
                             "clone_id",
                             {
                              data_type => "integer",
                              default_value => undef,
                              is_foreign_key => 1,
                              is_nullable => 0,
                              size => 4,
                             },
                            );
    __PACKAGE__->set_primary_key("clone_feature_id");
    __PACKAGE__->add_unique_constraint("clone_feature_clone_id_key", ["clone_id",'feature_id']);
    __PACKAGE__->belongs_to(
                            'feature',
                            'Bio::Chado::Schema::Sequence::Feature',
                            {
                             'foreign.feature_id' => 'self.feature_id' },
                           );

}

__END__

=head1 NAME

bsub_mirror_potato_bacs.pl - script to do something

=head1 SYNOPSIS

  bsub_mirror_potato_bacs.pl [options] args

  Options:

    -R <dir>
     set directory where the potato bac repository is kept
     will create directory if it does not exist.
     default /data/prod/public/potato_genome/bacs

    -v verbose
     print some information about what is being done

    -U update every BAC in the published file repository (to add
       metadata for example)

=head1 MAINTAINER

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 AUTHOR

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
