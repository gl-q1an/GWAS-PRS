# 1.1_qc_het.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript dir/1.1_qc_het.R inputfile outputfile")
} else {
  input <- args[1]
  output <- args[2]
}

het <- read.table(input, head=TRUE)
het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
write.table(het_fail[,1:2], output, row.names=FALSE, col.names=FALSE, quote=FALSE)