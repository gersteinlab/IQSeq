#!/bin/bash

gene=$1;
known_iso_only=KNOWN_ISO_ONLY;
iso=known_iso;
time_tag=`date +%Y%m%d-%H%M`
for s_cost in `seq 0 1 20`;
do
	sp_cost=`expr 20 - $s_cost`;
	./bin/sim_iso mle 0.01 200 $known_iso_only ~/bioinfo/data/hg18/$gene-isoforms.txt 0 $s_cost $sp_cost >> ~/bioinfo/data/isoform_mle/sim_iso/$gene-$iso-mle-$time_tag.rst;
done
