#!/bin/bash

# Here we got the merged bfile from `2_pps.sh` pps_merge_${bfilename}_hg37_1kg_indsnp.{bed,bim,fam} 
#  and the population file ppinfo.txt
# We further clustered finer groups (e.g., from EAS to CHB JPT CHS CDX KHV) 
#  based on the results of the previous population stratification step.

set -e

# --- set the path ---
work_dir=~/your/work/path
bfilename=bfilename
plink_dir=~/path/plink
gcta_dir=~/path/gcta
hg37_1kg_dir=~/Ref/hg37_1kg_phase3
script_dir=$(dirname "$0")

pop_select=("OWN" "CHB" "CHS" "CDX" "JPT" "KHV")

cd ${work_dir}/2_pps

# Create the population file to select for eligible individuals
awk 'NR>1 {print "0\t"$1,$6}' ${hg37_1kg_dir}/all_phase3.psam > ppinfo2.txt

> ppinfo_selected.txt
while read -r line; do
    pop_column=$(echo "$line" | awk '{print $3}')
    for pop in "${pop_select[@]}"; do
        if [ "$pop_column" == "$pop" ]; then
            echo "$line" >> ppinfo_selected.txt
            break
        fi
    done
done < ppinfo2.txt

awk '{print $1,$2"\tOWN"}' ${bfilename}_qc_done.fam >> ppinfo_selected.txt

awk '{print $1,$2}' ppinfo_selected.txt > ppinfo_selected_keep.txt

${plink_dir}/plink --bfile pps_merge_${bfilename}_hg37_1kg_indsnp \
--keep ppinfo_selected_keep.txt \
--make-bed \
--out pps_merge_${bfilename}_hg37_1kg_indsnp_selected

# Use the selected data to do the PCA
${gcta_dir}/gcta --bfile pps_merge_${bfilename}_hg37_1kg_indsnp_selected \
--make-grm \
--thread-num 5 \
--out tmp_grm

${gcta_dir}/gcta --grm tmp_grm \
--pca 20 \
--thread-num 5 \
--out pps_merge_${bfilename}_hg37_1kg_selected_pca20

awk '{print $2, $3, $4}' pps_merge_${bfilename}_hg37_1kg_selected_pca20.eigenvec > ${bfilename}_selected_pp_pcforplot.txt
awk '{print $2,$3}' ppinfo_selected.txt > ppinfo_selected_r.txt

Rscript ${script_dir}/2.1_pps_pca_plot.R ${bfilename}_selected_pp_pcforplot.txt ppinfo_selected_r.txt ${bfilename}_pps2.png

# The Populaion Stratification is done!
cat << GLHF
 /\_/\  ==========================
( o.o ) ||  2.PPS Extra is done! ||
 > ^ < / =========================
GLHF