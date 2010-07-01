package CXGN::TomatoGenome::CmdLine::Command::db_load;
use Moose;
use MooseX::Aliases;
use namespace::autoclean;
use Carp;

extends 'CXGN::TomatoGenome::CmdLine::Command';
   with 'CXGN::TomatoGenome::CmdLine::DBConnector';

sub abstract { <<'' }
reinitialize the genomic.clone_feature table, which links genomic.clone and the chado feature table



sub execute {
    my ( $self, $opt, $args ) = @_;

#TODO: truncate the clone_feature table and reload it with the most
#recent version of each bac seq that is present in the feature table

# old sql that did this without regard to sequence versions:
#  begin;

# select count(*) from genomic.clone_feature;

# truncate genomic.clone_feature;

# insert into genomic.clone_feature (clone_id,feature_id)
# select c.clone_id,f.feature_id
# from
#   (
#        select
# 	 f.feature_id as feature_id
#         ,f.name as name,
# 	,( case	when substring(f.name from 4 for 3) = 'HBa' then 'LE_HBa'
# 		when substring(f.name from 4 for 3) = 'SLe' then 'SL_EcoRI'
# 		when substring(f.name from 4 for 3) = 'SLm' then 'SL_MboI'
# 		end
# 	 ) as shortname
#         ,cast(substring(f.name from 7 for 4) as int) as platenum
# 	,substring(f.name from 11 for 1) as wellrow
# 	,cast(substring(f.name from 12 for 2) as smallint) as wellcol
#      from feature f
#      join cvterm c on (f.type_id = c.cvterm_id)
#      where f.name like 'C%.%' and c.name = 'BAC_clone'
#   ) as f
# join genomic.clone c
#   using(platenum,wellrow,wellcol)
# join genomic.library l
#   on(f.shortname = l.shortname and c.library_id = l.library_id)
# ;

# select count(*) from genomic.clone_feature;

# rollback;

}
