#!/bin/bash

# Here we got the plink bfile from `3_2impute.sh` ${bfilename}_impqcdone.{bed,bim,fam} 
#  and covar file from `2.3_pca.sh` ${bfilename}_covar.txt and covar file includes such as sex and age
#  you should merge the two files to `covar_file`
#  and a file includes FIDs IID and phenotypes

set -e

# --- set the path ---
work_dir=~/your/work/path
bfilename=bfile
plink_dir=~/path/plink
pheno_file=
covar_file=
ana_dir=${work_dir}/4_analysis

mkdir -p ${ana_dir}
cd ${ana_dir}

#####  Association Study #####
# If binary traits, use --logisitic --beta
# If continous traits, use --linear
# You can use --mpheno or --covar-name or --covar-number if necessary

${plink_dir}/plink --bfile ${work_dir}/3_imp/${bfilename}_impqcdone \
--pheno ${pheno_file} \
--logistic hide-covar \
--covar ${bfilename}_allcovar.txt \
--adjust \
--out result_${bfilename}

# Check the result
awk 'NR==1 || $12 < 5e-8' result_${bfilename}.assoc.linear > result_${bfilename}_sig.assoc.linear
samplenum=$(awk 'NR==2 {print $6}' result_${bfilename}.assoc.linear)
line=$(wc -l < result_${bfilename}_sig.assoc.linear)
sig=$((line - 1))
echo "The numble of sample to analysis is ${samplenum}" > result_${bfilename}.txt
echo "The numble of significant SNPs (p<5e-8) is ${sig}" > result_${bfilename}.txt

cat << GLHF
 /\_/\  ========================
( o.o ) ||  4 Analyse is done! ||
 > ^ < / =======================
GLHF