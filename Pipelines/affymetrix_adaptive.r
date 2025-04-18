# il codice è in grado di discriminare la tipologia di piattaforma specifica della tipologia Affymex utilizzata 
# sulla base della piattaforma, separa correttamente i dataset

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("oligo", quietly = TRUE))
  BiocManager::install("oligo")

library(oligo)
library(affy)

source("utils/geo_fetch_mapping.r")

# ----- funzione di utility per caricare l’annotation package corretto -------
load_platform_db <- function(platform) {
    # Determine appropriate annotation DB package
    if (grepl("^pd\\.", platform)) {
        # oligo annotation mapping: e.g. pd.hugene.1.0.st.v1 -> hugene10sttranscriptcluster.db
        parts <- strsplit(sub("^pd\\.", "", platform), "\\.")[[1]]
        organism <- parts[1]
        major <- parts[2]
        minor <- parts[3]
        # construct base name
        pkg_base <- paste0(organism, major, minor, "sttranscriptcluster")
        db_pkg <- paste0(pkg_base, ".db")
    } else {
        # affy annotation packages use <platform>.db
        db_pkg <- paste0(platform, ".db")
    }
    # install if missing
    if (!requireNamespace(db_pkg, quietly = TRUE)) {
        BiocManager::install(db_pkg, ask = FALSE, update = FALSE)
    }
    suppressPackageStartupMessages(library(db_pkg, character.only = TRUE))
    # annotation object name matches the DB package itself
    annot_db <- get(db_pkg)
    return(annot_db)
}



# Lanciare lo script da dentro ASD-BioMarkers
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

   
    cel_files <- list.files(absolute_path,
                            pattern = "\\.CEL(\\.gz)?$",
                            full.names = TRUE, recursive = TRUE)
  
    
    if (length(cel_files) == 0)
      stop("Nessun file CEL trovato in ", absolute_path)
    
    cat(length(cel_files), "CEL file(s) trovati.\n")

    # Determina i chip type
    cat("Analisi del tipo di chip per ciascun file...\n")
    oligo_compatible <- c()
    affy_compatible <- c()

    for (f in cel_files) {
      tryCatch({
        chip_type <- affyio::read.celfile.header(f)$cdfName
        if (!is.null(chip_type) && length(chip_type) == 1 && !is.na(chip_type) && nzchar(chip_type)) {
          if (grepl("HuGene|HTA|ST", chip_type, ignore.case = TRUE)) {
            cat(chip_type, "compatibile con oligo\n")
            oligo_compatible <- c(oligo_compatible, f)
          } else {
            cat(chip_type, "compatibile con affy\n")
            affy_compatible <- c(affy_compatible, f)
          }
          cat(basename(f), "->", chip_type, "\n")
        } else {
          cat("Chip type mancante o invalido per", basename(f), "\n")
        }
      }, error = function(e) {
        print(e)
        # cat("Errore nella lettura chip type:", f, "\n")
      })
    }

    cat("Check piattaforma eseguito\n")

    # Processa con oligo
    if (length(oligo_compatible) > 0) {
      cat("Normalizzazione con oligo...\n")
      celSet <- oligo::read.celfiles(oligo_compatible)
      # solo summarizzazione: niente background correction né normalizzazione
      eset <- oligo::rma(
          celSet,
          background = FALSE,
          normalize  = FALSE,
          target     = "core"
      )
      exprs_mat <- exprs(eset)

      probe_ids <- rownames(exprs_mat)


      platform   <- annotation(eset)       

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
                              paste0(basename(cel_dir_path), "_gene_expression_oligo.csv"))
      write.csv(agg_df, file = output_file, row.names = FALSE)
    
    }

    # Processa con affy
    if (length(affy_compatible) > 0) {
      cat("Normalizzazione con affy...\n")
      data_affy <- ReadAffy(filenames = affy_compatible)
      eset <- affy::rma(
          data_affy,
          background = FALSE,
          normalize  = FALSE,
          target     = "core"
      )
      exprs_mat <- exprs(eset)

      probe_ids <- rownames(exprs_mat)


      platform   <- annotation(eset)       

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
                              paste0(basename(cel_dir_path), "_gene_expression_affyo.csv"))
      write.csv(agg_df, file = output_file, row.names = FALSE)

    }

    cat("Conversione completata.\n") # nolint

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
