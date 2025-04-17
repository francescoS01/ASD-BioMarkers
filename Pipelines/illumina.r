 
 # Ensure BiocManager is available for Bioconductor packages
 if (!requireNamespace("BiocManager", quietly = TRUE)) {
     install.packages("BiocManager")
 }
 
 # List of Bioconductor packages needed
 bioc_pkgs <- c(
   "lumi",
   "lumiHumanAll.db",
   "lumiHumanIDMapping",
   "AnnotationDbi",
   "GEOquery",
   "Biobase",
   "DBI",
   "RSQLite"
 )
 
 # Install missing Bioconductor packages
 for (pkg in bioc_pkgs) {
     if (!requireNamespace(pkg, quietly = TRUE)) {
         BiocManager::install(pkg, ask = FALSE, update = FALSE)
     }
 }
 
 # List of CRAN packages needed
 cran_pkgs <- c("tools")
 
 # Install missing CRAN packages
 for (pkg in cran_pkgs) {
     if (!requireNamespace(pkg, quietly = TRUE)) {
         install.packages(pkg)
     }
 }
 
 # Load libraries
 library(BiocManager)
 library(lumi)
 library(lumiHumanAll.db)
 library(lumiHumanIDMapping)
 library(AnnotationDbi)
 library(GEOquery)
 library(Biobase)
 library(DBI)
 library(RSQLite)
 library(tools)
 


pipeline <- function(cel_dir_path) {
    cat("Starting Illumina pipeline for folder:", cel_dir_path, "\n")
    
    root <- file.path("Datasets", "Illumina")

    absolute_path = file.path(root, cel_dir_path)

    parent_dir <- dirname(absolute_path)
    if (!dir.exists(parent_dir)) {
        dir.create(parent_dir, recursive = TRUE)
    }

    # Create CSV directory inside the GEO folder
    csv_dir_path <- file.path(absolute_path, "CSV")
    if (!dir.exists(csv_dir_path)) {
        dir.create(csv_dir_path, recursive = TRUE)
    }

    cat(absolute_path,"\n")
   
    files <- list.files(absolute_path, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
    cat(length(files), "\n")
    
    for (input_file in files) {
      cat("Processo file ",input_file,"\n")
      # derive a base name and output path for each file
      file_base <- tools::file_path_sans_ext(basename(input_file))
      output_file <- file.path(csv_dir_path, paste0(file_base, "_gene_expression.csv"))

      # Read Illumina data
      # lumiData <- lumiR(input_file,)

      lumiData <- tryCatch({
        lumiR(input_file, probeID = "ID_REF", sep = "\t", header = TRUE)
      }, error = function(e) {
        cat("Error reading file:", input_file, "\n")
        cat("Error message:", e$message, "\n")
        return(NULL) # Return NULL if there's an error reading the file
      })

      

      # Annotate LumiBatch with nuID and gene symbol
      lumiData <- addNuID2lumi(lumiData, lib.mapping="lumiHumanIDMapping")

      # Extract expression matrix
      exprs_mat <- exprs(lumiData)

      probe_ids <- rownames(exprs_mat)
      annot <- AnnotationDbi::select(lumiHumanAll.db,
                               keys = probe_ids,
                               columns = c("PROBEID", "SYMBOL"),
                               keytype = "PROBEID")

      # Merge annotation (SYMBOL) with expression matrix
      exprs_df <- as.data.frame(exprs_mat)
      exprs_df$PROBEID <- probe_ids
      exprs_df$SYMBOL <- annot$SYMBOL[match(probe_ids, annot$PROBEID)]
      merged <- exprs_df

      agg_df <- aggregate(
        . ~ SYMBOL,
        data = merged[, !(names(merged) %in% "PROBEID")],
        FUN = mean,
        na.rm = TRUE
    )
      # Export to CSV: each row is a probe with its gene symbol and intensities
      write.csv(agg_df,
              file = output_file,
              row.names = FALSE)
    }
}

# Execute pipeline based on command-line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
    cat("Running pipeline on dataset:", args[1], "\n")
    pipeline(args[1])

} else {
    stop("Usage: Rscript illumina.r <dataset_folder_name>\n")
}
