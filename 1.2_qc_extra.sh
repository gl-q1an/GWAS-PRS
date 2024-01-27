#!/bin/bash

# Here we got the plink bfile from `1_qc.sh` ${bfilename}_qc2_done.{bed,bim,fam}

set -e

# --- set the path ---
work_dir=~/your/work/path
qc_dir=${work_dir}/1_qc
bfilename=bfilename
plink_dir=~/path/plink
plink2_dir=~/path/plink2
annotref_file=~/Ref/chrpostors/snp150_hg19.txt

cd ${qc_dir}

##### STEP Extra1 Remove Duplicated #####
${plink2_dir}/plink2 --bfile ${bfilename}_qc2_done \
--set-all-var-ids @:# \
--make-bed \
--out ${bfilename}_qc3_1duptmp1

${plink2_dir}/plink2 --bfile ${bfilename}_qc3_1duptmp1 \
--rm-dup force-first \
--make-bed \
--out ${bfilename}_qc3_1nodup

rm ${bfilename}_qc3_1duptmp1.{bed,bim,fam}
##### STEP Extra1 is DONE! q(>_<)9 #####

##### STEP Extra2 Annotation #####
${plink_dir}/plink --bfile ${bfilename}_qc3_1nodup \
--update-name ${annotref_file} \
--make-bed \
--out ${bfilename}_qc3_2annot
##### STEP Extra2 is DONE! q(>_<)9 #####

# The QC Extra is done!
cat << GLHF
 /\_/\  ========================
( o.o ) ||  1.QC Extrais done! ||
 > ^ < / =======================
GLHF