# GWAS and PRS
Here is a tutorial on how to conduct Genome-wide Association Study (GWAS). I added some details to make the whole process more holistic.
In addition, how to use genotyping data from UKB for GWAS analysis is included in this tutorials. See the `README_UKB.md` for detail. The following tutorials are for genotyping data.

# Tools
---
There are many tools to choose from for PCA, Phase and Impute, it is said that the speed and results may not be exactly the same from one tool to another, but I haven't actually done it in different ways so I can only give examples of the ones I've used, **_feel free to give me feedback if you have a better tool to recommend_!**

- plink1.9/plink2 (for qc and association): https://www.cog-genomics.org/plink
- KING (for famliy relation): https://www.kingrelatedness.com/
- GTCA (for population stratification): https://yanglab.westlake.edu.cn/software/gcta/
- EIGENSOFT (smartpca for PCA): https://www.hsph.harvard.edu/alkes-price/software/
- Eagle (for phase): https://alkesgroup.broadinstitute.org/Eagle/
- Minimac4 (for impute): https://genome.sph.umich.edu/wiki/Minimac4
- bgzip & tabix: https://www.htslib.org/download/
- R 4.2.3 (for data orgnization and visualization)

# Data Prepare
---
I don't know which kind of data format you got, there may be multiple batches that need to be merged, and there may be data in different formats that need to be converted, but all in all we need to convert it to a set of plink1.9 binary format (`.bed`, `.bim`, `.fam`) for QC. Familiarity with the operation of plink as well as R programming may be helpful in this step, here is an example of our own data, unfortunately it can not be uploaded, if you want to practice you can try the data in the tutorial(https://github.com/MareesAT/GWA_tutorial).

The data I got is plink binary data, no data conversion needed.

```bash
# See how many SNPs and individuals
wc -l yourfile.bim yourfile.fam
# See the files
head yourfile.fam
head yourfile.bim
tail yourfile.bim
head -n 50000 yourfile.bim | tail
```

## Check how many SNPs and individuals
When examining the number of SNPs and individuals, the number of SNPs for chr 1-22 sequenced by microarray are generally around 700,000 (my file is 743,722) and you need to pay attention to whether the number of individuals matches the number of individuals you sent for testing.
## Check the fam files
When examining the fam file, you need to focus on FID, IID, gender, and phenotype. Here is an example from my fam file.  Here the IID number is the chip location number.

```
0 20461111111_R11C11 0 0 2 -9
```

The file does not have FID, and there is a "\_" in the IID (I remember the time when I processed the data before, some software would take the part in front of the "\_" as the FID, and after the "\_" as the IID, which led to the error in my later analysis), so we'd better change the FID to the same, for example, the name of the place + the serial number "wh001". 

```bash
awk '{print $1"\t"$2"\twh"sprintf("%04d",NR)"\twh"sprintf("%04d",NR)}' your.fam > rawID2proID.txt
plink --bfile yourfile --update-ids rawID2proID.txt --make-bed --out pro
```

Sex is not missing, if it is missing it may need to be added based on general new documentation or sex chromosomes, or this data someone else has already done sex QC and it doesn't affect your next step in the analysis.
The phenotype is missing, so in the following QC (e.g., hwe), we can't QC according to case/control separately.
## Check the bim files
When examining the bim file, you need to check which chromosomes the bim file includes and how the loci are named.

```
# head yourfile.bim
0	10:100301704	0	0	0	G
...
# tail yourfile.bim
26	ilmnseq_rs386829316	0	16527	T	C
...
# head -n 50000 yourfile.bim | tail
1	rs10494997	0	215460047	A	G
1	ilmnseq_1:215460832	0	215460832	T	C
...
```

We can find that bim chromosomes range from 0-26, 0 may refer to SNPs that are missing some of the information (I'm not sure), 26 refers to mitochondrial chromosomes, and we only need data for chromosomes 1-22 in most cases. In addition, we can see that SNPs on chromosomes 1-22 are named differently, and we need to ANNOTE later if needed.
## Data merging
In fact, it's possible to have your sequencing company merge it for you, which I did do. If you merge on your own, you may run into chain-flip problems, which you can solve using `plink --flip`, as described in https://www.cog-genomics.org/plink2/data#merge3. If you have multiple files to merge, it is recommended to generate a `.missnp` file chain of issues for each file before merging with `plink --merge-list`.

# Quality Control
---
According to the idea of ​​this article "Data quality control in genetic case-control association studies"(http://www.ncbi.nlm.nih.gov/pubmed/21085122), perform individual quality control first and then perform SNPs quality control to preserve SNPs as much as possible.  Personally, I don't think the order has much effect, as you can see, the article mentioned at the beginning was published after this one, but not strictly in this order.
- individual level
	- missing rate > 0.03
	- sex mismatch
	- heterozygosity deviated ± 3 stand deviation (s.d.)
	- family related (by KING, 3 degree)
- SNPs level
	- call rate < 0.97
	- minor allele frequency (MAF) < 0.01
	- Hardy-Weinberg equilibrium (P < 1e-10)
	- NOT biallelic SNPs
In some references, the criterion for HW equilibrium is P < 1e-6 in controls and P < 1e-10 in cases. However, considering that at the very beginning we may not have written the phenotypic information into the fam file and did not necessarily do a dichotomous phenotype, P<1e-10 was used for all.
## QC extra

Sometimes our data is not that ideal and additional QC may be needed after the above QC is done, such as removal of duplicate sites, annotation to rsID (these steps may not be needed if your data is very ideal).

When we browsed the bim file before, we found that SNPs loci are named in different forms, such as rsID, chip, kgpID, etc., which may have duplicate IDs.

>Note: The reason we use plink2 to convert SNPID to chr:pos and then remove the duplications is because the data we use to annotate doesn't have A1 information, and we have to make sure that there are no duplicate sites on the position information. plink's `--list-duplicate-vars suppress-first`  based on the position and A1, and it can't recognize positional duplicates but A1 is not duplicated.

After removing duplicates, here we conduct the annotation using the file `snp150_hg19.txt` like: 

```
chromosome:start        name
1:10039 rs978760828
1:10043 rs1008829651
...
```

# Population Stratification & PCA
---
You may wonder why I use two different software when both steps are essentially PCA, because my teacher taught me that `GCTA` is faster, so it's easier to calculate with this method after merging with 1KG data, while `engensoft` is more accurate in calculating, so it's more suitable for principal component calculation. In addition to these two tools, `plink` can also perform PCA.

When performing these steps, we need to remove the areas of high LD, referring to the website https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)

>Note: When I was under Ubuntu 18.06, R 3.6+ and engensoft 7 were conflicting and couldn't be installed at the same time, so I had to install engensoft under the virtual environment conda.

## Populaion Stratification

As mentioned earlier, population stratification is the process of merging your own data with data from various ethnic groups from 1000 Genomic Project and then running PCA to see what demographic information your data clusters with. Data stratified by different populations need to be analyzed separately and may need to be rejected if the sample size is small. If you're sure of the ethnicity of your sample, you can skip this step and go straight to PCA (e.g. I'm sure all my data is from the Chinese Han population).

Our gene structure is hg37, and we download the 1000 Genomic Project Data from https://www.cog-genomics.org/plink/2.0/resources. 

By checking the information in the files, We found that the population was divided into 5 major categories (AFR, AMR, EAS, EUR, SAS) and 26 sub-categories, We can do PCA on the larger population first and then pull out the smaller populations within it as the tutorial did. Our data is aggregated with East Asian populations and we further select East Asian populations for population stratification. Or, you can just based on the needs of your data, e.g., if most of our data is expected to be CHB data, we combine the CEU, FIN, GBR, etc. populations together as EUR, and the CHB, CHS, and JPT are separate. 

>Important: Make sure your variants are named by rsID, as the 1000G data is using the rsID and the two need to be merged later. If you can't find `snp150_hg19.txt`, you can use the `pvar` file or the `bim` file after QC from 1000G to make a file similar to `snp150_hg19.txt` for annotation.

## PCA

There are many different ways to select the number of principal components, one of which can be through the p-value of Tracy-Widom statistics, for example here we select the first 5 principal components. This step  generates a file of covariates and excludes individuals outside of 6 standard deviations at the same time.

# Phase & Impute

The template and processing of eagle can be referred to the official website, which also mentions that if your sample is twice as large as the template ,the reference is not necessary.

Minimac4 reference can be downloaded from the official website.

After Imputation, similar to the SNPs QC above, quality control and annotation of SNPs is still required. INFO represents the quality of impute, generally choose INFO>0.8
# Association Analysis

If it is a binary variable use `--logistic`, if it is a continuous variable use `--linear`, if there are covariates you can add `hide-covar` after, otherwise  each covariate will be calculated separately, and the final generated file will be very large.

Manhattan plot without going into detail, you can use CMplot (https://github.com/YinLiLin/CMplot), you can also refer to ggplot2 Manhattan drawing tutorials. The main difficulty in plotting using ggplot2 is how to determine the x-axis based on chromosome length, which is addressed in this tutorial(https://r-graph-gallery.com/101_Manhattan_plot.html), `Manhattan` package can not adjust the color so it is not beautiful enough.
# PRS

...
# DATA from UK Biobank

The phenotype data in UK Biobank has already been partially processed, so it is different from the way we process our own data. See the `README_UKB.md` for detail.

>The reason why I write this airticle is that the genotype of data of my lab needs to be processed.  So until the ukb data needs to be processed, I may not write this part. (I'm so lazy.) So if you have any questions, You can contact me to discuss it, after all I'm also a learner. 
