begin;

truncate clone_feature;

insert into clone_feature (clone_id,feature_id) select c.clone_id,f.feature_id
from
(select
	f.feature_id as feature_id,
	f.name as name,
	( case	when substring(f.name from 4 for 3) = 'HBa' then 'LE_HBa'
		when substring(f.name from 4 for 3) = 'SLe' then 'SL_EcoRI'
		when substring(f.name from 4 for 3) = 'SLm' then 'SL_MboI'
		end
	) as shortname,
	cast(substring(f.name from 7 for 4) as int) as platenum,
	substring(f.name from 11 for 1) as wellrow,
	cast(substring(f.name from 12 for 2) as smallint) as wellcol
from feature f
join cvterm c on (f.type_id = c.cvterm_id)
where f.name like 'C%.%' and c.name = 'BAC_clone'
) as f
join genomic_bt.clone c
  using(platenum,wellrow,wellcol)
join genomic_bt.library l
  on(f.shortname = l.shortname and c.library_id = l.library_id)
;

commit;
