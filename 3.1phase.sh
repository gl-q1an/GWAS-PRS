#!/bin/bash

# Here we got the plink bfile from `2.3_pca.sh` ${bfilename}_ppspca_done.{bed,bim,fam} 
#  and ${bfilename}_indep_snplist.prune.in
# And Get 1000 Genomic file from https://www.cog-genomics.org/plink/2.0/resources

set -e

# --- set the path ---
work_dir=~/your/work/path
bfilename=bfilename
plink_dir=~/path/plink
eagle_dir=~/path/eigensoft/bin/
bgzip_dir=~/path/htslib/bin
tabix_dir=~/path/htslib/bin
samtools_dir=~/path/samtools/bin
bcftools_dir=~/path/bcftools/bin
eagle_ref=~/Ref/eagleref/hg37_v5a_eagle
eagle_ref_raw=~/Ref/1kg_ref/hg37_phase3_v5a

phase_dir=${work_dir}/3_phase
mkdir -p ${phase_dir}
cd ${phase_dir}
cp ${work_dir}/2_pps/${bfilename}_ppspca_done.{bed,bim,fam} ${phase_dir}/

##### STEP1 Phase #####
for chr in {1..22}; do
${plink_dir}/plink --bfile ${bfilename}_ppspca_done \
--chr ${chr} \
--recode vcf-iid \
--out ${bfilename}_forphase_chr${chr}

${bgzip_dir}/bgzip ${bfilename}_forphase_chr${chr}.vcf
${tabix_dir}/tabix -p vcf ${bfilename}_forphase_chr${chr}.vcf.gz
done

# ====== Generate reference, if necessary =====
# Files in ${eagle_ref_raw} are as follows:
# human_g1k_v37.fasta
# ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
# ALL.chr{1..22}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
#  genetic_map_hg19_withX.txt.gz from https://alkesgroup.broadinstitute.org/Eagle/downloads/tables/
cd ${eagle_ref}
cp ${eagle_ref_raw}/human_g1k_v37.fasta ./

${samtools_dir}/samtools faidx human_g1k_v37.fasta

for chr in {1..22}; do
  ${bcftools_dir}/bcftools view --no-version -Ou -c 2 ${eagle_ref_raw}/ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | \
  ${bcftools_dir}/bcftools norm --no-version -Ou -m -any | \
  ${bcftools_dir}/bcftools norm --no-version -Ob -o ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf -d none -f human_g1k_v37.fasta && \
  ${bcftools_dir}/bcftools index -f ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf
done
# Warning: '-Ob' means compressed bcf file, but the file is named .bcf.
# Because eagle can only identified '.bcf' rather than '.bcf.gz' file.
# Although there seems to be no difference in the results
#==============================================

cd ${phase_dir}
for chr in {1..22};do
${eagle_dir}/eagle \
--vcfRef ${eagle_ref}/ALL.chr${chr}.phase3_integrated.20130502.genotypes.bcf \
--vcfTarget ${bfilename}_forphase_chr${chr}.vcf.gz \
--geneticMapFile ${eagle_ref}/genetic_map_hg19_withX.txt.gz \
--numThreads 15 \
--vcfOutFormat=z \
--outPrefix ${bfilename}_phased_chr${chr}
done
##### STEP 1 is DONE! q(>_<)9 #####

# The Phase is done!
cat << GLHF
 /\_/\  ========================
( o.o ) ||  3 1Phase is done! ||
 > ^ < / =======================
GLHF