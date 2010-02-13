set search_path=sgn_people,public;

--------------------------------
-- CLONE IL MAPPING PROJECT LOG
--------------------------------

drop trigger sp_project_il_mapping_clone_log_uniqueness on sp_project_il_mapping_clone_log;
drop table sp_project_il_mapping_clone_log;

CREATE TABLE sp_project_il_mapping_clone_log (
    sp_project_il_mapping_clone_log_id serial primary key NOT NULL,
    sp_project_id integer references sp_project,
    sp_person_id integer references sp_person,
    clone_id integer references genomic_bt.clone,
    is_current boolean DEFAULT true,
    created timestamp without time zone DEFAULT now()
   
);

comment on table  sp_project_il_mapping_clone_log is 'linking table showing which sp_project is currently assigned to map a given clone to the zamir IL lines.  also provides a modification history with is_current and created columns';


CREATE OR REPLACE FUNCTION check_sp_project_il_mapping_clone_log_uniqueness ()
        RETURNS trigger AS $$
	$_[0] = 'sgn_people.sp_project_il_mapping_clone_log';
	$_[1] = 'clone_id';

        if($_TD->{new}{is_current} ne 'f'){

                # make sure there are no other current
                my $r = spi_exec_query("SELECT * FROM $_[0] WHERE $_[1] = $_TD->{new}{$_[1]} and is_current = true");

                if ($r->{processed} > 0){
                        die "Can't insert a second current assignment of $_[1] $_TD->{new}{$_[1]}";
                }
        }

        # everything ok so far
        return;

$$ LANGUAGE plperl;

CREATE TRIGGER sp_project_il_mapping_clone_log_uniqueness
        BEFORE INSERT
        ON sgn_people.sp_project_il_mapping_clone_log
        FOR EACH ROW
        EXECUTE PROCEDURE check_sp_project_il_mapping_clone_log_uniqueness();


grant select on  sgn_people.sp_project_il_mapping_clone_log to public;
grant insert on  sgn_people.sp_project_il_mapping_clone_log to web_usr;
grant update on  sgn_people.sp_project_il_mapping_clone_log to web_usr;
grant select on  sgn_people.sp_project_il_mapping_clone_l_sp_project_il_mapping_clone_l_seq to web_usr;
grant update on  sgn_people.sp_project_il_mapping_clone_l_sp_project_il_mapping_clone_l_seq to web_usr;

--------------------------------
-- CLONE IL MAPPING BIN LOG
--------------------------------

drop trigger clone_il_mapping_bin_log_uniqueness on sp_clone_il_mapping_bin_log;
drop table clone_il_mapping_bin_log;

CREATE TABLE clone_il_mapping_bin_log (
	sp_clone_il_mapping_bin_log_id serial primary key NOT NULL,
	genotype_region_id integer references phenome.genotype_region,
	sp_person_id integer references sgn_people.sp_person,
	clone_id integer references genomic_bt.clone,
	is_current boolean DEFAULT true,
	created timestamp without time zone DEFAULT now()
);

comment on table sgn_people.clone_il_mapping_bin_log is 'linking table showing which phenome.genotype_region a given clone has been mapped to. also provides a modification history with its is_current and created columns';

CREATE OR REPLACE FUNCTION check_clone_il_mapping_bin_log_uniqueness()
        RETURNS trigger AS $$
	$_[0] = 'sgn_people.clone_il_mapping_bin_log';
	$_[1] = 'clone_id';

        if($_TD->{new}{is_current} ne 'f'){

                # make sure there are no other current
                my $r = spi_exec_query("SELECT * FROM $_[0] WHERE $_[1] = $_TD->{new}{$_[1]} and is_current = true");

                if ($r->{processed} > 0){
                        die "Can't insert a second current assignment of $_[1] $_TD->{new}{$_[1]}";
                }
        }

        # everything ok so far
        return;

$$ LANGUAGE plperl;

CREATE TRIGGER clone_il_mapping_bin_log_uniqueness
        BEFORE INSERT
        ON sgn_people.clone_il_mapping_bin_log
        FOR EACH ROW
        EXECUTE PROCEDURE check_clone_il_mapping_bin_log_uniqueness();

grant select on  clone_il_mapping_bin_log to public;
grant insert on  clone_il_mapping_bin_log to web_usr;
grant update on  clone_il_mapping_bin_log to web_usr;
grant select on  clone_il_mapping_bin_log_sp_clone_il_mapping_bin_log_id_seq to web_usr;
grant update on  clone_il_mapping_bin_log_sp_clone_il_mapping_bin_log_id_seq to web_usr;


--------------------------------
-- CLONE VALIDATION LOG
--------------------------------

drop trigger clone_verification_log_uniqueness on clone_verification_log;
drop table clone_verification_log;

CREATE TABLE clone_verification_log (
    clone_verification_log_id serial NOT NULL,
    sp_person_id integer references sp_person,
    clone_id integer references genomic_bt.clone,
    ver_int_read boolean NOT NULL DEFAULT false,
    ver_bac_end boolean NOT NULL DEFAULT false,
    is_current boolean DEFAULT true,
    created timestamp without time zone DEFAULT now()
);

comment on table clone_verification_log is 'table showing which clones have been validated by a variety of methods.  columns may be added to this without warning.  details about each validation experiment should be written into the comment field on the detail page for the clone';

CREATE OR REPLACE FUNCTION check_clone_verification_log_uniqueness()
        RETURNS trigger AS $$
	$_[0] = 'sgn_people.clone_verification_log';
	$_[1] = 'clone_id';

        if($_TD->{new}{is_current} ne 'f'){

                # make sure there are no other current
                my $r = spi_exec_query("SELECT * FROM $_[0] WHERE $_[1] = $_TD->{new}{$_[1]} and is_current = true");

                if ($r->{processed} > 0){
                        die "Can't insert a second current assignment of $_[1] $_TD->{new}{$_[1]}";
                }
        }

        # everything ok so far
        return;

$$ LANGUAGE plperl;

CREATE TRIGGER sp_clone_verification_log_uniqueness
        BEFORE INSERT
        ON sgn_people.clone_verification_log
        FOR EACH ROW
        EXECUTE PROCEDURE check_clone_verification_log_uniqueness();

grant select on  sgn_people.clone_verification_log to public;
grant insert on  sgn_people.clone_verification_log to web_usr;
grant update on  sgn_people.clone_verification_log to web_usr;
grant select on  sgn_people.clone_verification_log_clone_verification_log_id_seq to web_usr;
grant update on  sgn_people.clone_verification_log_clone_verification_log_id_seq to web_usr;


