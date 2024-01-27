# 2.1_pps_pca_plot.R

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  vecfile <- args[1]
  output <- args[2]
} else if (length(args) == 3) {
  vecfile <- args[1]
  popfile <- args[2]
  output <- args[3]
} else {
  stop("Usage: Rscript dir/2.1_pps_pca_plot.R File_with_ID_PC1_PC2 population.info(Optional) outputffigure")
}

if (!require("ggplot2", character.only = TRUE)) {
  stop("ggplot2 is necessary!")
} else {
  suppressMessages(library(ggplot2))
}

vec <- read.table(vecfile, header=FALSE)
colnames(vec) <- c("ID","PC1","PC2")

if (exists("popfile")) {
  pop <- read.table(popfile, header=FALSE)
  colnames(pop) <- c("ID","POP")
  vec_pop <- merge(vec, pop, by="ID")
  
  p <- ggplot(vec_pop, aes(x = PC1, y = PC2, color = POP)) +
    geom_point(size=2,alpha=0.8,stroke=0) +
    labs(title = "Populaion Stratification Plot", x = "PC1", y = "PC2") +
    theme_bw() +
    theme( 
      legend.position = "right",
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5),
    )
} else {
  p <- ggplot(vec, aes(x = PC1, y = PC2, color='red')) +
    geom_point(size=2,alpha=0.8,stroke=0) +
    labs(title = "PCA Plot", x = "PC1", y = "PC2") +
    theme_bw() +
    theme(
      legend.position="none", 
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5),
    )
}

ggsave(output, plot = p, width = 10, height = 7)
cat(paste0("Congratulations, ",output," is finished!!!"))