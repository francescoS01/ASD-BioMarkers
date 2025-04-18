# Ensure BiocManager is installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")

if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery", ask = FALSE)
suppressPackageStartupMessages(library(GEOquery))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))

pipeline <- function(agilent_dir) {
    cat("Starting Agilent pipeline for folder:", agilent_dir, "\n")
    csv_dir <- file.path(agilent_dir, "pipelineOut")
    if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)
    files <- list.files(agilent_dir, pattern = "\\.txt$", full.names = TRUE)
    expr_list <- list()
    for (file in files) {
        cat("Processing file:", file, "\n")
        # load file skipping metadata
        lines <- readLines(file)
        hdr <- grep("^FEATURES", lines)[1]
        if (is.na(hdr)) stop("Cannot find FEATURES header in: ", file)
        df <- read.delim(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = hdr - 1, check.names = FALSE)
        if (!all(c("SystematicName", "gProcessedSignal") %in% colnames(df))) {
            stop("Required columns missing in file: ", file)
        }
        sample_code <- tools::file_path_sans_ext(basename(file))
        expr <- aggregate(gProcessedSignal ~ SystematicName,
                          data = df,
                          FUN = mean,
                          na.rm = TRUE)
        names(expr)[2] <- sample_code
        expr_list[[1 + length(expr_list)]] <- expr
    }
    # merge all samples by SystematicName using base merge
    if (length(expr_list) == 0) stop("No data frames to merge")
    final_mat <- expr_list[[1]]
    if (length(expr_list) > 1) {
        for (i in 2:length(expr_list)) {
            final_mat <- merge(final_mat, expr_list[[i]], by = "SystematicName", all = TRUE)
        }
    }
    # set ProbeID from original probe names
    final_mat$ProbeID <- as.character(final_mat$SystematicName)
    final_mat$SystematicName <- NULL

    # retrieve platform annotation for GPL6480
    gpl <- getGEO("GPL6480", AnnotGPL=TRUE)
    anno_df <- Table(gpl)
    cols <- colnames(anno_df)
    id_col <- cols[1]
    sym_candidates <- grep("Symbol", cols, ignore.case = TRUE, value = TRUE)
    if (length(sym_candidates) == 0) stop("Gene symbol column not found in GPL6480 annotation")
    sym_col <- sym_candidates[1]
    anno <- anno_df[, c(id_col, sym_col)]
    colnames(anno) <- c("ProbeID", "Symbol")
    anno$ProbeID <- as.character(anno$ProbeID)

    # join and aggregate by Gene Symbol
    df_sym <- merge(final_mat, anno, by = "ProbeID", all.x = TRUE)
    # fill missing gene symbols with probe ID
    df_sym$Symbol[is.na(df_sym$Symbol)] <- df_sym$ProbeID
    expr_sym <- aggregate(. ~ Symbol, data = df_sym[ , !(names(df_sym) %in% "ProbeID")], FUN = mean, na.rm = TRUE)

    # write expression matrix with Gene Symbol
    write.csv(expr_sym,
              file = file.path(csv_dir, "expression_matrix.csv"),
              row.names = FALSE)
}

# run from command line
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript agilent.r <dataset_code>")
}
codename <- args[1]
root <- file.path("Datasets", "Agilent")
agilent_dir <- file.path(root, codename)

pipeline(agilent_dir)