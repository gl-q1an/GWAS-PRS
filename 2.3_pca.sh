#!/bin/bash

# Here we got the plink bfile from `1_pps.sh` ${bfilename}_qc_done.{bed,bim,fam} 
#  and ${bfilename}_indep_snplist.prune.in
# And Get 1000 Genomic file from https://www.cog-genomics.org/plink/2.0/resources

set -e

# --- set the path ---
work_dir=~/your/work/path
bfilename=bfilename
plink_dir=~/path/plink
smartpca_dir=~/path/eigensoft/bin/
script_dir=$(dirname "$0")

cd ${work_dir}/2_pps

##### STEP 1 Writing and generating parameter related files #####
${plink_dir}/plink --bfile ${bfilename}_qc_done \
--extract ${bfilename}_indep_snplist.prune.in \
--recode \
--out ${bfilename}_indep_pca

echo "##transfer.conf
genotypename: ${bfilename}_indep_pca.ped
snpname: ${bfilename}_indep_pca.map
indivname: ${bfilename}_indep_pca.ped
outputformat: EIGENSTART
genotypeoutname: ${bfilename}_indep_pca.geno
snpoutname: ${bfilename}_indep_pca.snp
indivoutname: ${bfilename}_indep_pca.ind
familynames: NO" >${bfilename}_transfer.conf

echo "##runningpca.conf
genotypename: ${bfilename}_indep_pca.geno
snpname: ${bfilename}_indep_pca.snp
indivname: ${bfilename}_indep_pca.ind
evecoutname: ${bfilename}_indep_pca.evec
evaloutname: ${bfilename}_indep_pca.eval
altnormstyle: NO
numoutevec: 20
numoutlieriter: 5
outliersigmathresh: 6.0" >${bfilename}_runningpca.conf
##### STEP 1 is DONE! q(>_<)9 #####

##### STEP 2 Running Pca and plot #####
${smartpca_dir}/convertf -p ${bfilename}_transfer.conf
${smartpca_dir}/smartpca -p ${bfilename}_runningpca.conf > ${bfilename}_smartpca.log

awk 'NR>1 {print $1,$2,$3}' ${bfilename}_indep_pca.evec > ${bfilename}_smartpca_rplot.txt

Rscript ${script_dir}/2.1_pps_pca.R ${bfilename}_smartpca_rplot.txt ${bfilename}_pca.png
##### STEP 2 is DONE! q(>_<)9 #####

##### STEP 3 Generate covariate files and remove individuals beyond 6 s.d. #####
awk 'NR>1 {print $1,$2,$3,$4,$5,$6}' ${bfilename}_indep_pca.evec > ${bfilename}_covar.txt
awk 'NR>1 {print $1,$1}' ${bfilename}_indep_pca.evec > ${bfilename}_pca_ind_keep.txt

${plink_dir}/plink --bfile ${bfilename}_qc_done \
--keep ${bfilename}_pca_ind_keep.txt \
--make-bed \
--out ${bfilename}_ppspca_done
##### STEP 3 is DONE! q(>_<)9 #####

# The PCA is done!
cat << GLHF
 /\_/\  ===========================
( o.o ) ||  2.PPS & PCA is done! ||
 > ^ < / ==========================
GLHF