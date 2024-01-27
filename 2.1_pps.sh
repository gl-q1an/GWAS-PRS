#!/bin/bash

# Here we got the plink bfile from `1_qc.sh` ${bfilename}_qc2_done.{bed,bim,fam} 
#  or `1.2_qc_extra.sh` ${bfilename}_qc3_2annot.{bed,bim,fam}
# And Get 1000 Genomic file from https://www.cog-genomics.org/plink/2.0/resources

set -e

# --- set the path ---
work_dir=~/your/work/path
bfilename=bfilename
plink_dir=~/path/plink
plink2_dir=~/path/plink2
gcta_dir=~/path/gcta
hg37_1kg_dir=~/Ref/hg37_1kg_phase3
highLD_hg37_file=~/Ref/GWAS/highLD/37highLD-regions.txt
script_dir=$(dirname "$0")

pps_dir=${work_dir}/2_pps
mkdir -p ${pps_dir}
# choose your file  ${bfilename}_qc2_done.{bed,bim,fam} or ${bfilename}_qc3_2annot.{bed,bim,fam}
bfileprefix=${bfilename}_qc3_2annot
cp ${work_dir}/1_qc/${bfileprefix}.bed ${pps_dir}/${bfilename}_qc_done.bed
cp ${work_dir}/1_qc/${bfileprefix}.bim ${pps_dir}/${bfilename}_qc_done.bim
cp ${work_dir}/1_qc/${bfileprefix}.fam ${pps_dir}/${bfilename}_qc_done.fam

##### STEP 1 Prepare the data from 1000 Genimic #####
# If you have already processed this data, you can proceed directly to the next step.
#===== 1.1 Download the 1000 Genomic Data =====
# https://www.cog-genomics.org/plink/2.0/resources
mkdir -p ${hg37_1kg_dir}
cd ${hg37_1kg_dir}
# all_phase3.pgen.zst(2.25GB)
wget https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst
# all_phase3.pvar.zst(1.26GB)
wget https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst
# phase3_corrected.psam
wget https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam
# deg2_phase3.king.cutoff.out.id
wget https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id
# De-compress all_phase3.pgen(6.3G) all_phase3.pvar(12G)
${plink2_dir}/plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
${plink2_dir}/plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar

#==== 1.2 Prepare for the data from 1kg =====
mv phase3_corrected.psam all_phase3.psam

# Check the population information in the files
# awk 'NR>1 {print $5}' all_phase3.psam | sort | uniq -c
# awk 'NR>1 {print $6}' all_phase3.psam | sort | uniq -c

${plink2_dir}/plink2 --pfile all_phase3 \
--remove deg2_phase3.king.cutoff.out.id \
--rm-dup exclude-mismatch \
--set-missing-var-ids '@:#:$1:$2' \
--chr 1-22 \
--mind 0.03 \
--geno 0.03 \
--maf 0.05 \
--snps-only just-acgt \
--max-alleles 2 \
--make-bed \
--out hg37_1kg_qc
##### STEP 1 is DONE! q(>_<)9 #####

##### STEP 2 Prepare the data from 1000 Genimic and ours #####
#===== 2.1 Extract the same SNPs =====
cd ${pps_dir}

${plink_dir}/plink --bfile ${bfilename}_qc_done \
--indep-pairwise 50 5 0.2 \
--exclude ${highLD_hg37_file} \
--range \
--out ${bfilename}_indep_snplist

${plink_dir}/plink --bfile ${hg37_1kg_dir}/hg37_1kg_qc \
--extract ${bfilename}_indep_snplist.prune.in \
--make-bed \
--out pps_hg37_1kg_indsnp

awk '{print$2}' pps_hg37_1kg_indsnp.bim > pps_${bfilename}_hg37_1kg.snplist

${plink_dir}/plink --bfile ${bfilename}_qc_done \
--extract pps_${bfilename}_hg37_1kg.snplist \
--make-bed \
--out pps_${bfilename}_indsnp

#===== 2.2 Check the strand and merge =====
${plink_dir}/plink --bfile pps_${bfilename}_indsnp \
--bmerge pps_hg37_1kg_indsnp.bed pps_hg37_1kg_indsnp.bim pps_hg37_1kg_indsnp.fam \
--merge-mode 6 \
--out pps_1mismatch_${bfilename}_hg37_1kg

if [ -s "pps_1mismatch_${bfilename}_hg37_1kg.missnp" ]; then
# if `pps_1mismatch_${bfilename}_hg37_1kg.missnp` is not empty
    ${plink_dir}/plink --bfile pps_hg37_1kg_indsnp \
    --flip pps_1mismatch_${bfilename}_hg37_1kg.missnp \
    --make-bed \
    --out pps_hg37_1kg_indsnp_flip

    ${plink_dir}/plink --bfile pps_${bfilename}_indsnp \
    --bmerge pps_hg37_1kg_indsnp_flip.bed pps_hg37_1kg_indsnp_flip.bim pps_hg37_1kg_indsnp_flip.fam \
    --merge-mode 6 \
    --out pps_2mismatch_${bfilename}_hg37_1kg_flip1kg

    # exclude variants that do not complete the chain conversion
    ${plink_dir}/plink --bfile pps_${bfilename}_indsnp \
    --exclude pps_2mismatch_wh2023_hg37_1kg_flip1kg.missnp \
    --make-bed \
    --out pps_${bfilename}_indsnp_merge

    ${plink_dir}/plink --bfile pps_hg37_1kg_indsnp \
    --exclude pps_2mismatch_wh2023_hg37_1kg_flip1kg.missnp \
    --make-bed \
    --out pps_hg37_1kg_indsnp_merge

    # merge
    ${plink_dir}/plink --bfile pps_${bfilename}_indsnp_merge \
    --bmerge pps_hg37_1kg_indsnp_merge.bed pps_hg37_1kg_indsnp_merge.bim pps_hg37_1kg_indsnp_merge.fam \
    --make-bed \
    --out pps_merge_${bfilename}_hg37_1kg_indsnp
else
    # If there's no issue with strand, just merge!
    ${plink_dir}/plink --bfile pps_${bfilename}_indsnp \
    --bmerge pps_hg37_1kg_indsnp.bed pps_hg37_1kg_indsnp.bim pps_hg37_1kg_indsnp.fam \
    --make-bed \
    --out pps_merge_${bfilename}_hg37_1kg_indsnp
fi
##### STEP 2 is DONE! q(>_<)9 #####

##### STEP 3 Population Stafication #####
# Use of data from all populations
${gcta_dir}/gcta --bfile pps_merge_${bfilename}_hg37_1kg_indsnp \
--make-grm \
--thread-num 5 \
--out tmp_grm

${gcta_dir}/gcta --grm tmp_grm \
--pca 20 \
--thread-num 5 \
--out pps_merge_${bfilename}_hg37_1kg_pca20

awk '{print $2, $3, $4}' pps_merge_${bfilename}_hg37_1kg_pca20.eigenvec > ${bfilename}_allpp_pcforplot.txt

awk 'NR>1 {print $1,$5}' ${hg37_1kg_dir}/all_phase3.psam > ppinfo.txt
awk '{print $2"\tOWN"}' ${bfilename}_qc_done.fam >> ppinfo.txt

Rscript ${script_dir}/2.1_pps_pca_plot.R ${bfilename}_allpp_pcforplot.txt ppinfo.txt ${bfilename}_pps1.png

##### STEP 3 is DONE! q(>_<)9 #####
# Using data from a subset of the population.
# See `2.2_pps_extra.sh`

# The QC is done!
cat << GLHF
 /\_/\  ====================
( o.o ) ||  2.PPS is done! ||
 > ^ < / ===================
GLHF