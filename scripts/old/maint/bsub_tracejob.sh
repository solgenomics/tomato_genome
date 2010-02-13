#!/bin/sh
echo ----- looking up the job in the processing log -----
echo select timestamp,host,message from cxgn_bac_pipeline_processing_log where message like "'%$1%'" order by timestamp desc limit 2 | psql -U web_usr -h db.sgn.cornell.edu cxgn
log=`ls /data/shared/tomato_genome/bsub_annotate_all/job_logs/$1.*`
echo ----- catting job log $log ------------------------------
sudo cat $log
echo ----- running tracejob, finding execution node ---
tracejob $1 | tee /tmp/bsub_tracejob.tmp
host= `perl -e 'while(<>) { if(/(\w+.cluster.sgn)/) { print "$1\n"; last }}' /tmp/bsub_tracejob.tmp`
echo ----- sshing to execution host $host, catting job record -
ssh $host cat /var/spool/torque/spool/$1*
echo ----- sshing to execution host $host, grepping mom logs --
sudo ssh $host grep $1 /var/spool/torque/mom_logs/*



