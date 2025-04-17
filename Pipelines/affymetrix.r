if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("oligo", quietly = TRUE))
  BiocManager::install("oligo")

  
library(oligo)