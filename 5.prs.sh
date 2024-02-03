#!/bin/bash

# Here we got the plink bfile from `3_2impute.sh` ${bfilename}_impqcdone.{bed,bim,fam} 
#  and covar file from `2.3_pca.sh` ${bfilename}_covar.txt and covar file includes such as sex and age
#  you should merge the two files to `covar_file` to calculate the PRS
# Base summary file is also needed!

set -e 

# --- set the config ---
work_dir=
script_dir=$(dirname "$0")
bfilename=
plink_dir=
prsice_dir=
prscs_dir=
prscs_ref=
base_summary=${work_dir}/0_rawdata/base_summary.txt.gz
base_pheno=
pheno_file=
covar_file=

n_gwas=

prs_dir=${work_dir}/5_prs
mkdir -p ${prs_dir}
cd ${prs_dir}

python ${script_dir}/format_summary.py \
 --summary ${base_summary} \
 --out ${base_pheno}_std.txt.gz \
 --outform SNP,CHR,BP,A1,A2,Frq,P,BETA,SE,N,INFO

##### STEP1. QC of the base summary and target #####
# If necessary, perform build conversions by liftover.
gunzip -c ${base_pheno}_std.txt.gz |\
awk 'NR==1 || ($6 > 0.01 && $6 < 0.99) {print}' |\
awk 'NR==1 || ($4 ~ /^(A|T|C|G)$/ && $5 ~ /^(A|T|C|G)$/) {print}' |\
gzip > ${base_pheno}_std_1snpqc.txt.gz

# This step may take miniutes
gunzip -c ${base_pheno}_std_1snpqc.txt.gz |\
awk '{seen[$1]++; if(seen[$1]==1){print}}' |\
gzip > ${base_pheno}_std_2nodup.txt.gz

# Palindromic SNPs
gunzip -c ${base_pheno}_std_2nodup.txt.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > ${base_pheno}_std_3baseqc.txt.gz

# Mismatching SNPs
# It is usually not required!
##### STEP 1 is DONE! q(>_<)9 #####

##### STEP2. Calculate the PRS by PRSice2 #####
# The clump will be done defaultly!
# 1.Calculate the best P value, and the ${bfilename}_${base_pheno}_prs.best has the PRS best fit
Rscript ${prsice_dir}/PRSice.R \
 --prsice ${prsice_dir}/PRSice_linux \
 --base ${base_pheno}_std_3baseqc.txt.gz \
 --target ${work_dir}/3_imp/${bfilename}_impqcdone \
 --clump-kb 250kb \
 --clump-r2 0.1 \
 --clump-p 1 \
 --pheno ${pheno_file} \
 --pheno-col state \
 --cov ${covar_file} \
 --cov-col sex,age,@PC[1-5] \
 --binary-target T \
 --beta \
 --thread 15 \
 --out ${prs_dir}/${bfilename}_${base_pheno}_prs

# 2. Calculate the PRS for all samples
Rscript ${prsice_dir}/PRSice.R \
 --prsice ${prsice_dir}/PRSice_linux \
 --base ${base_pheno}_std_3baseqc.txt.gz \
 --target ${work_dir}/3_imp/${bfilename}_impqcdone \
 --clump-kb 250kb \
 --clump-r2 0.1 \
 --clump-p 1 \
 --no-regress \
 --bar-levels 0.01,0.05,0.1 \
 --fastscore \
 --beta \
 --thread 15 \
 --out ${prs_dir}/${bfilename}_${base_pheno}_prs_samples
##### STEP 2 is DONE! q(>_<)9 #####

##### STEP3. Calculate the PRS by PRS-cs-auto #####
# Calculating PRS using different methods can be used for sensitivity analysis.
gunzip -c ${base_pheno}_std_3baseqc.txt.gz | awk '{print $1"\t"$4"\t"$5"\t"$8"\t"$9}' > ${base_pheno}_prscs_sum.txt

python ${prscs_dir}/PRScs.py \
 --ref_dir=${prscs_ref} \
 --bim_prefix=${work_dir}/3_imp/${bfilename}_impqcdone \
 --sst_file=${base_pheno}_prscs_sum.txt \
 --n_gwas=${n_gwas} \
 --out_dir=${prs_dir}/${bfilename}_${base_pheno}_prscs

mkdir -p ${prs_dir}/${bfilename}_${base_pheno}_prscs_chr
> ${prs_dir}/${bfilename}_${base_pheno}_prscs_pst_eff_a1_b0.5_phiauto_all.txt
for i in {1..22};do
cat ${prs_dir}/${bfilename}_${base_pheno}_prscs_pst_eff_a1_b0.5_phiauto_chr${i}.txt >> ${prs_dir}/${bfilename}_${base_pheno}_prscs_pst_eff_a1_b0.5_phiauto_all.txt
mv ${prs_dir}/${bfilename}_${base_pheno}_prscs_pst_eff_a1_b0.5_phiauto_chr${i}.txt ${prs_dir}/${bfilename}_${base_pheno}_prscs_chr/
done

${plink_dir}/plink --bfile ${work_dir}/3_imp/${bfilename}_impqcdone \
 --score ${prs_dir}/${bfilename}_${base_pheno}_prscs_pst_eff_a1_b0.5_phiauto_all.txt 2 4 6 \
 --out ${prs_dir}/${bfilename}_${base_pheno}_prscs_res
##### STEP 3 is DONE! q(>_<)9 #####

cat << GLHF
 /\_/\  ========================
( o.o ) ||   5 PRS is done!   ||
 > ^ < / =======================
GLHF