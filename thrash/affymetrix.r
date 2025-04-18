if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

bioc_pkgs <- c(
    "affy",            # pacchetto classico 3'/U133
    "oligo",           # per Gene/Transcript‑ST compatibili con oligo
    "affyio",          # per leggere header CEL
    "AnnotationDbi",   # mapping probe‑>gene
    "GEOquery",        # metadati GEO
    "Biobase", "DBI", "RSQLite"
)

for (pkg in bioc_pkgs)
    if (!requireNamespace(pkg, quietly = TRUE))
        BiocManager::install(pkg, ask = FALSE, update = FALSE)

cran_pkgs <- c("tools")
for (pkg in cran_pkgs)
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg)

library(affy)
library(AnnotationDbi)
library(GEOquery)
library(Biobase)
library(DBI)
library(RSQLite)
library(tools)
library(affyio)   # accesso rapido a header CEL
library(oligo)    # processamento Gene/Transcript‑ST

# ----- funzione di utility per caricare l’annotation package corretto -------
load_platform_db <- function(platform) {
    # es. "hgu133plus2" -> "hgu133plus2.db"
    db_pkg <- paste0(platform, ".db")
    if (!requireNamespace(db_pkg, quietly = TRUE))
        BiocManager::install(db_pkg, ask = FALSE, update = FALSE)
    suppressPackageStartupMessages(library(db_pkg, character.only = TRUE))
    get(db_pkg)   # restituisce l’oggetto AnnotationDb
}

# ----- eventuale funzione esterna per mapping campioni‑>diagnosi ------------
source("utils/geo_fetch_mapping.r")   # già usata in illumina_v2.r

# ----------------------------- pipeline -------------------------------------
pipeline <- function(cel_dir_path) {

    cat("Starting Affymetrix pipeline for folder:", cel_dir_path, "\n")
    root          <- file.path("Datasets", "Affymetrix")
    absolute_path <- file.path(root, cel_dir_path)

    if (!dir.exists(absolute_path))
        stop("Directory non trovata: ", absolute_path)

    # cartella di output
    csv_dir_path <- file.path(absolute_path, "pipelineOut")
    if (!dir.exists(csv_dir_path))
        dir.create(csv_dir_path, recursive = TRUE)

    # -------- lettura dei CEL ------------------------------------------------
    cel_files <- list.files(absolute_path,
                            pattern = "\\.CEL(\\.gz)?$",
                            full.names = TRUE, recursive = TRUE)
    
    if (length(cel_files) == 0)
        stop("Nessun file CEL trovato in ", absolute_path)
    
    cat(length(cel_files), "CEL file(s) trovati.\n")
    
    # ---- determina chip‑type per scegliere affy vs oligo --------------------
    chip_type <- affyio::read.celfile.header(cel_files[1])$cdfName
    is_oligo  <- grepl("HuGene|HTA|ST", chip_type, ignore.case = TRUE)
    
    if (is_oligo) {
        cat("Chip riconosciuto come compatibile con 'oligo' (", chip_type, ").\n")
        celSet <- oligo::read.celfiles(cel_files)
        # solo summarizzazione: niente background correction né normalizzazione
        eset <- oligo::rma(
            celSet,
            background = FALSE,
            normalize  = FALSE,
            target     = "core"
        )
    } else {
        cat("Chip riconosciuto come 3'/U133‑style (", chip_type, "). Uso 'affy'.\n")
        affyData <- ReadAffy(filenames = cel_files)
        eset <- expresso(
            affyData,
            bg.correct       = FALSE,
            normalize        = FALSE,
            pmcorrect.method = "pmonly",
            summary.method   = "avgdiff"
        )
    }
    
    exprs_mat <- exprs(eset)

    probe_ids <- rownames(exprs_mat)


    platform   <- annotation(eset)       
    cat("Piattaforma riconosciuta:", platform, "\n")
    annot_db   <- load_platform_db(platform)

    annot <- AnnotationDbi::select(
        annot_db,
        keys     = probe_ids,
        columns  = c("PROBEID", "SYMBOL"),
        keytype  = "PROBEID"
    )

 
    exprs_df         <- as.data.frame(exprs_mat)
    exprs_df$PROBEID <- probe_ids
    exprs_df$SYMBOL  <- annot$SYMBOL[match(probe_ids, annot$PROBEID)]

    merged <- exprs_df

    agg_df <- aggregate(
        . ~ SYMBOL,
        data = merged[, !(names(merged) %in% "PROBEID")],
        FUN  = mean,          # cambia qui se vuoi max/median…
        na.rm = TRUE
    )

    # -------- export ---------------------------------------------------------
    output_file <- file.path(csv_dir_path,
                             paste0(basename(cel_dir_path), "_gene_expression.csv"))
    write.csv(agg_df, file = output_file, row.names = FALSE)
    cat("Gene‑level expression salvata in:", output_file, "\n")


    cat("Running geo_fetch_mapping …\n")
    geo_output_path <- csv_dir_path
    get_Info_from_GEO(cel_dir_path, geo_output_path)
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 1) {
    pipeline(args[1])
} else {
    stop("Uso: Rscript affymetrix.r <folder_dataset>\n")
}