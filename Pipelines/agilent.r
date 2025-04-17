if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("lumi", quietly = TRUE))
  BiocManager::install("lumi")

  
library(limma)
