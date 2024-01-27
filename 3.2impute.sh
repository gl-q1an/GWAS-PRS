#!/bin/bash

# Here we got the plink bfile from `3_1phase.sh` ${bfilename}_chr${chr}_phased.vcf.gz 

set -e

# --- set the path ---
work_dir=~/your/work/path
phase_dir=${work_dir}/3_phase
bfilename=bfile
plink_dir=~/path/plink
plink2_dir=~/path/plink2
minimac4_dir=~/path/Minimac4/release-build
minimac4_ref=~/Ref/minimac4ref
annotref_file=~/Ref/chrpostors/snp150_hg19.txt

geno_thresh=0.03
maf_thresh=0.01
hwe_thresh=1e-10

imp_dir=${work_dir}/3_imp
mkdir -p ${imp_dir}
cd ${imp_dir}

##### STEP 1 Impute #####
for chr in {1..22};do
${minimac4_dir}/minimac4 \
--refHap ${minimac4_ref}/${chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
--haps ${phase_dir}/${bfilename}_phased_chr${chr}.vcf.gz \
--prefix ${bfilename}_chr${chr}_imputed \
--cpus 15
done
##### STEP 1 is DONE! q(>_<)9 #####

##### STEP 2 QC after Imputation #####
#===== INFO > 0.8 (R-squared Rsq > 0.8) =====
for chr in {1..22};do
awk 'NR>1 && $7 > 0.8 {print $1,$2,$7}' ${bfilename}_chr${chr}_imputed.info > ${bfilename}_chr${chr}_imputed_info08.txt
awk '{print $1}' ${bfilename}_chr${chr}_imputed_info08.txt > ${bfilename}_imputed_chr${chr}_info08_snplist.txt
echo "chr${chr} is done!"
done

#===== vcf to bed =====
for chr in {1..22};do
${plink_dir}/plink --vcf ${bfilename}_chr${chr}_imputed.dose.vcf.gz \
--extract ${bfilename}_imputed_chr${chr}_info08_snplist.txt \
--allow-no-sex \
--double-id \
--make-bed \
--out ${bfilename}_imputed_chr${chr}_info08
done

#===== Combination =====
> ${bfilename}_mergelist.txt
for chr in {2..22};do
echo ${bfilename}_imputed_chr${chr}_info08 >> ${bfilename}_mergelist.txt
done

${plink_dir}/plink --bfile ${bfilename}_imputed_chr1_info08 \
--merge-list ${bfilename}_mergelist.txt \
--hwe ${hwe_thresh} \
--maf ${maf_thresh} \
--geno ${geno_thresh} \
--snps-only just-acgt \
--allow-no-sex \
--make-bed \
--out ${bfilename}_impdone

#===== Annotation =====
${plink2_dir}/plink2 --bfile ${bfilename}_impdone \
--set-all-var-ids @:# \
--make-bed \
--out ${bfilename}_impdone_tmp

${plink2_dir}/plink2 --bfile ${bfilename}_impdone_tmp \
--rm-dup force-first \
--make-bed \
--out ${bfilename}_impdone_tmp_nodup

rm -f ${bfilename}_impdone_tmp.{bed,bim,fam}

${plink_dir}/plink --bfile ${bfilename}_impdone_tmp_nodup \
--update-name ${annotref_file} \
--make-bed \
--out ${bfilename}_impqcdone

rm -f ${bfilename}_impdone_tmp_nodup.{bed,bim,fam}
##### STEP2 is DONE! q(>_<)9 #####

# The Impute is done!
cat << GLHF
 /\_/\  ========================
( o.o ) ||  3 Impute is done! ||
 > ^ < / =======================
GLHF