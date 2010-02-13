set search_path=genomic;

select count(*) as "Total BAC-End chromatograms" from chromat;

select count(*) as "Total basecalled bac end sequences" from gss;

select library.shortname as "Library", count(*) as "BACs with 2 ends >= 150 hqi and no problems"
from clone as cn 
join library
  using(library_id)
join chromat as c1 
  ON (c1.clone_id=cn.clone_id AND c1.primer='T7') 
join gss as g1 
  ON ( c1.chromat_id = g1.chromat_id ) 
join chromat as c2 
  ON(cn.clone_id=c2.clone_id AND c2.primer='SP6') 
join gss as g2 
  ON ( c2.chromat_id = g2.chromat_id ) 
WHERE c1.chromat_id IS NOT NULL 
  and c2.chromat_id IS NOT NULL 
  and g2.flags = 0 
  and g1.flags = 0
group by library.shortname;

select library.shortname as "Library", count(*) as "BACs with only T7 >= 150 hqi and no problems" 
from library
join clone as cn
  using(library_id)
join chromat as c1 
  ON (c1.clone_id=cn.clone_id AND c1.primer='T7') 
join gss as g1 
  ON ( c1.chromat_id = g1.chromat_id ) 
left join chromat as c2 
  ON(cn.clone_id=c2.clone_id AND c2.primer='SP6') 
left join gss as g2 
  ON ( c2.chromat_id = g2.chromat_id ) 
WHERE g1.flags = 0
  and (
    c2.chromat_id IS NULL 
    OR g2.flags != 0 
  )
group by library.shortname;

select library.shortname as "Library",count(*) as "BACs with only SP6 >= 150 hqi and no problems"
from library
join clone as cn 
  using(library_id)
left join chromat as c1 
  ON (c1.clone_id=cn.clone_id AND c1.primer='T7') 
left join gss as g1 
  ON ( c1.chromat_id = g1.chromat_id ) 
join chromat as c2 
  ON(cn.clone_id=c2.clone_id AND c2.primer='SP6') 
join gss as g2 
  ON ( c2.chromat_id = g2.chromat_id ) 
WHERE (c1.chromat_id IS NULL 
    OR g1.flags = 0
  )	
  and g2.flags = 0 
  and c2.chromat_id IS NOT NULL 
group by library.shortname;

select count(*) as "Total BAC ends with >= 150 hqi and no problems"
from gss as g 
where g.flags=0;

----------------------------------------------------------
--------------- contamination ----------------------------
----------------------------------------------------------

select library.shortname as "Library", count(*) as "BACs with two contaminated ends"
from library
join clone as cn
  using(library_id)
join chromat as c1
  ON (c1.clone_id=cn.clone_id AND c1.primer='T7')
join gss as g1
  ON ( c1.chromat_id = g1.chromat_id )
join chromat as c2
  ON(cn.clone_id=c2.clone_id AND c2.primer='SP6')
join gss as g2
  ON ( c2.chromat_id = g2.chromat_id )
WHERE g1.flags&32 !=0 AND g2.flags&32 != 0
group by library.shortname;


select library.shortname as "Library", count(*) as "BACs with only T7 end showing contamination"
from library
join clone as cn
  using(library_id)
join chromat as c1
  ON (c1.clone_id=cn.clone_id AND c1.primer='T7')
join gss as g1
  ON ( c1.chromat_id = g1.chromat_id )
join chromat as c2
  ON(cn.clone_id=c2.clone_id AND c2.primer='SP6')
join gss as g2
  ON ( c2.chromat_id = g2.chromat_id )
WHERE g1.flags&32 !=0 AND g2.flags&32 = 0
group by library.shortname;

select library.shortname as "Library", count(*) as "BACs with only SP6 end showing contamination"
from clone as cn
join library
  using(library_id)
join chromat as c1
  ON (c1.clone_id=cn.clone_id AND c1.primer='T7')
join gss as g1
  ON ( c1.chromat_id = g1.chromat_id )
join chromat as c2
  ON(cn.clone_id=c2.clone_id AND c2.primer='SP6')
join gss as g2
  ON ( c2.chromat_id = g2.chromat_id )
WHERE g1.flags&32 = 0 AND g2.flags&32 != 0
group by library.shortname;

select l.shortname as "Library", count(*) as "Total BAC ends showing contamination" 
from gss as g 
join qc_report as q
  using(gss_id)
join chromat as chr
  using(chromat_id)
join clone as cn
  using(clone_id)
join library as l
  using(library_id)
where q.hqi_length >= 100 
  and g.flags&32 != 0
group by l.shortname;

