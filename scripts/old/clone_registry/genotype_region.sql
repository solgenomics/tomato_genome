set search_path = phenome;

drop table genotype_region;
create table genotype_region (
	genotype_region_id serial primary key,
	genotype_id int references genotype,
	marker_id_nn bigint references sgn_bt.marker (marker_id),
	marker_id_ns bigint not null references sgn_bt.marker (marker_id),
	marker_id_sn bigint not null references sgn_bt.marker (marker_id),
	marker_id_ss bigint references sgn_bt.marker (marker_id),
	zygocity_code varchar(1) check( zygocity_code in ('a','b','c','d','h') ),
	lg_id int not null references sgn_bt.linkage_group,
	type varchar(32) not null check( type in ('bin','map','inbred') ),
	name varchar(32),
	sp_person_id int references sgn_people.sp_person,
	modified_date timestamp with time zone default now() not null,
	create_date timestamp with time zone default now() not null,
	obsolete boolean not null default false
);

comment on table genotype_region
is 'polymorphic regions from a genotype, delineated by markers in a certain linkage group on a certain map';
comment on column genotype_region.genotype_id
is 'optional genotype this region belongs to.  some regions are artificial, arising from combinations of other regions, and thus do not have an associated genotype';
comment on column genotype_region.lg_id
is 'the linkage group in a specific version of a specific map where this region is located';
comment on column genotype_region.type
is 'the type of polymorphic region this is.  map is mapping experiment data, inbred is IL lines segments, and bin is a derived region based on a boolean combination of inbred fragments.  For bin regions, the specific boolean combination of fragments that make the bin is not stored.';
comment on column genotype_region.sp_person_id
is 'the person who loaded this datum.  optional';
comment on column genotype_region.name
is 'special name for this region, if any.  optional';
comment on column genotype_region.marker_id_nn
is 'the north marker in the pair of markers bracketing the north end of this region. this may be null for regions at the north end of a linkage group';
comment on column genotype_region.marker_id_ns
is 'the south marker in the pair of markers bracketing the north end of this region';
comment on column genotype_region.marker_id_sn
is 'the north marker in the pair of markers bracketing the south end of this region';
comment on column genotype_region.marker_id_ss
is 'the south marker in the pair of markers bracketing the south end of this region. this may be null for regions at the south end of a linkage group.';


grant select on  phenome.genotype_region to public;
