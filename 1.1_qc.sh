#!/bin/bash

# This script is for quality control
# Good luck and have fun!

set -e

# --- set the path ---
work_dir=~/your/work/path
plink_bfile=~/set/your/bplinkfilepath
plink_dir=~/path/plink
king_dir=~/path/king
script_dir=$(dirname "$0")
bfilename=$(basename ${plink_bfile})

# --- set the parameter ---
mind_thresh=0.03
geno_thresh=0.03
maf_thresh=0.01
hwe_thresh=1e-10

# Create the qc directory
qc_dir=${work_dir}/1_qc
mkdir -p ${qc_dir}
cp ${plink_bfile} ${qc_dir}/
cd ${qc_dir} || exit

##### STEP 1 Individual level #####
#==== 1.1 Perform sex control and leave only autosomes =====
${plink_dir}/plink --bfile ${bfilename} \
--check-sex \
--out ${bfilename}_qc1_ind1_sex

grep "PROBLEM" ${bfilename}_qc1_ind1_sex.sexcheck |  awk '{print$1,$2}' > ${bfilename}_qc1_ind1_sexfail.txt

${plink_dir}/plink --bfile ${bfilename} \
--remove ${bfilename}_qc1_ind1_sexfail.txt \
--chr 1-22 \
--make-bed \
--out ${bfilename}_qc1_ind1_sexdone

#==== 1.2 Perform missing rate/het/family-related =====
#------- 1.2.1 missing rate --------
${plink_dir}/plink --bfile ${bfilename}_qc1_ind1_sexdone \
--missing \
--out ${bfilename}_qc1_ind2_mis

sed '1d' ${bfilename}_qc1_ind2_mis.imiss | awk -v thresh="${mind_thresh}" '$6 > thresh {print$1,$2}' > ${bfilename}_qc1_ind2_misfail.txt

#------- 1.2.2 het --------
${plink_dir}/plink --bfile ${bfilename}_qc1_ind1_sexdone \
--het \
--out ${bfilename}_qc1_ind3_het

Rscript ${script_dir}/1.1_qc_het.R ${bfilename}_qc1_ind3_het.het ${bfilename}_qc1_ind3_hetfail.txt

cat ${bfilename}_qc1_ind2_misfail.txt ${bfilename}_qc1_ind3_hetfail.txt | sort | uniq > ${bfilename}_qc1_ind3_mishetfail.txt

${plink_dir}/plink --bfile ${bfilename}_qc1_ind1_sexdone \
--remove ${bfilename}_qc1_ind3_mishetfail.txt \
--make-bed \
--out ${bfilename}_qc1_ind3_mishetdone

#------- 1.2.3 family related --------
${plink_dir}/plink --bfile ${bfilename}_qc1_ind3_mishetdone \
--indep-pairwise 50 5 0.2 \
--out ${bfilename}_qc1_ind4_indepSNP

${plink_dir}/plink --bfile ${bfilename}_qc1_ind3_mishetdone \
--extract ${bfilename}_qc1_ind4_indepSNP.prune.in \
--make-bed \
--out ${bfilename}_qc1_ind4_indep

${king_dir}/king -b ${bfilename}_qc1_ind4_indep.bed \
--unrelated --degree 3 --prefix ${bfilename}_qc1_ind4_un_

rm -f ${bfilename}_qc1_ind4_indep.{bim,fam,bed} ${bfilename}_qc1_ind4_indepSNP.prune*
# rm -f ${bfilename}_qc1_ind1_sexdone.{bim,fam,bed} ${bfilename}_qc1_ind3_mishetdone.{bim,fam,bed}

${plink_dir}/plink --bfile ${bfilename}_qc1_ind3_mishetdone \
--remove ${bfilename}_qc1_ind4_un_unrelated_toberemoved.txt \
--make-bed \
--out ${bfilename}_qc1_done
##### STEP 1 is DONE! q(>_<)9 #####

##### STEP 2 SNPs level #####
${plink_dir}/plink --bfile ${bfilename}_qc1_done \
--geno ${geno_thresh} \
--maf ${maf_thresh} \
--hwe ${hwe_thresh} \
--snps-only just-acgt \
--write-snplist \
--out ${bfilename}_qc2_snplist

${plink_dir}/plink --bfile ${bfilename}_qc1_done \
--extract ${bfilename}_qc2_snplist.snplist \
--make-bed \
--out ${bfilename}_qc2_done
##### STEP 2 is DONE! q(>_<)9 #####

# The QC is done!
cat << GLHF
 /\_/\  ====================
( o.o ) ||  1.QC is done! ||
 > ^ < / ===================
GLHF