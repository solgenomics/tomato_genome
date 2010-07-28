begin;

set search_path = genomic,public;

--\o /tmp/before.txt

--select
--      cf.clone_id
--     ,f.feature_id
--     , substring(f.name from 1 for 13) as name
--     , cast(substring(f.name from 15) as smallint) as v
--     , timelastmodified
--from clone_feature cf
--join feature f using(feature_id)
--join cvterm c on (f.type_id = c.cvterm_id)
--where
--   f.name like 'C%.%'
--   and f.name not like '%-%'
--order by
--    cf.clone_id,
--    f.timelastmodified desc,
--    cast(substring(f.name from 15) as smallint) desc -- highest version
--;

--\o

select
       distinct on(c.clone_id)
        c.clone_id
       ,f.feature_id
       , substring(f.name from 1 for 13) as name
       , cast(substring(f.name from 15) as smallint) as v
       , timelastmodified
into temporary table tomato_clone_feature
   from (
         select
   	f.feature_id as feature_id,
   	f.name as name,
           f.timelastmodified,
   	( case	when substring(f.name from 4 for 3) = 'HBa' then 'LE_HBa'
   		when substring(f.name from 4 for 3) = 'SLe' then 'SL_EcoRI'
   		when substring(f.name from 4 for 3) = 'SLm' then 'SL_MboI'
   		when substring(f.name from 4 for 3) = 'SLf' then 'SL_FOS'
   		end
   	) as shortname,
   	cast(substring(f.name from 7 for 4) as int) as platenum,
   	substring(f.name from 11 for 1) as wellrow,
   	cast(substring(f.name from 12 for 2) as smallint) as wellcol
         from feature f
         join cvterm c on (f.type_id = c.cvterm_id)
         where f.name like 'C%.%' and c.name = 'BAC_clone' and f.name not like '%-%'
        ) as f
   join clone c
     using(platenum,wellrow,wellcol)
   join library l
     on(f.shortname = l.shortname and c.library_id = l.library_id)
   order by
       c.clone_id,
       f.timelastmodified desc, -- most recent feature
       cast(substring(f.name from 15) as smallint) desc -- highest version
   ;

--\o /tmp/after.txt

--select * from tomato_clone_feature;

--\o

--select * 
--into temporary table clone_feature_backup
--from clone_feature;

delete from clone_feature
  where
       feature_id in ( select feature_id from tomato_clone_feature )
    or clone_id   in ( select clone_id   from tomato_clone_feature )
;

insert into clone_feature (feature_id,clone_id)
  select feature_id,clone_id from tomato_clone_feature;

--\o /tmp/deleted_clones.txt

select cf.clone_id,l.shortname || c.platenum || c.wellrow || c.wellcol
from clone_feature_backup cf
join clone c using(clone_id)
join library l using(library_id)
where cf.clone_id not in( select clone_id from clone_feature );

--\o

--rollback;
commit;
