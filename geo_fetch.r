if (!requireNamespace("GEOquery", quietly=TRUE))
    BiocManager::install("GEOquery")

library(GEOquery)
# Increase download timeout and set method for GEOquery
options(timeout = 600)
options(download.file.method.GEOquery = "libcurl")

get_Info_from_GEO <- function(geo){

  gse <- getGEO(geo, GSEMatrix=TRUE)
  eset <- gse[[1]]

  # estrai la metadata table
  pdata <- pData(eset)

  # Identify all 'characteristics' columns
  char_cols <- grep("^characteristics", colnames(pdata), value = TRUE)

  print(char_cols)
  # Extract diagnosis for each sample
  dx <- sapply(seq_len(nrow(pdata)), function(i) {
      # combine all characteristic strings for sample i
      char_vals <- unlist(strsplit(paste(pdata[i, char_cols], collapse = ";"), ";"))
      char_vals <- trimws(char_vals)
      # find any field starting with 'dx' or 'diagnosis'
      dx_entries <- char_vals[grepl("^(dx|diagnosis|disease)", char_vals, ignore.case = TRUE)]
      if (length(dx_entries) == 0) {
          return(NA_character_)
      } else {
          # remove everything up to the first colon and return the value
          sub("^[^:]+:\\s*", "", dx_entries[1])
      }
  })

  # costruiamo la tabella di mapping
  mapping <- data.frame(
    SampleID = rownames(pdata),
    Title     = as.character(pdata$title),
    Diagnosis = dx,
    stringsAsFactors = FALSE
  )

  # 6) dà un’occhiata
  head(mapping)
  # oppure salvala su disco
  
  write.csv(mapping, file = paste0(geo, "_diagnosis_mapping.csv"), row.names = FALSE)
}


# Execute pipeline based on command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
    cat("Running geo fetcher on dataset:", args[1], "\n")
    get_Info_from_GEO(args[1])

} else {
    stop("Usage: Rscript geo_fetch.r <geo>\n")
}
