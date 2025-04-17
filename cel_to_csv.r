# il codice è in grado di discriminare la tipologia di piattaforma specifica della tipologia Affymex utilizzata  # nolint
# sulla base della piattaforma, separa correttamente i dataset

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("oligo", quietly = TRUE))
  BiocManager::install("oligo")

library(oligo)
library(affy)

# Lanciare lo script da dentro ASD-BioMarkers
convert_cel_to_csv <- function(cel_dir_path) {

    root = "Datasets/Affymetrix_U133_Plus_v2"
    cat("Avvio della conversione...\n")
    
    absolute_path = file.path(root, cel_dir_path)
    # Crea il path alla cartella CSV
    parent_dir <- dirname(absolute_path)
    if (!dir.exists(parent_dir)) {
        dir.create(parent_dir, recursive = TRUE)
    }
    csv_dir_path <- file.path(parent_dir, "CSV")
    if (!dir.exists(csv_dir_path)) {
        dir.create(csv_dir_path)
    }
    
    # Lista di file .CEL
    cel_files <- list.celfiles(absolute_path, full.names = TRUE)

    
    if (length(cel_files) == 0) {
        cat("Nessun file .CEL trovato in:", absolute_path, "\n")
        return(NULL)
    }

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
      data_oligo <- read.celfiles(oligo_compatible)
      matrice_oligo <- exprs(oligo::rma(data_oligo))
      write.csv(matrice_oligo, file = file.path(csv_dir_path, "data_oligo.csv"))
    }

    # Processa con affy
    if (length(affy_compatible) > 0) {
      cat("Normalizzazione con affy...\n")
      data_affy <- ReadAffy(filenames = basename(affy_compatible), celfile.path = absolute_path)
      matrice_affy <- exprs(affy::rma(data_affy))
      write.csv(matrice_affy, file = file.path(csv_dir_path, "data_affy.csv"))
    }

    cat("Conversione completata.\n") # nolint
}

convert_cel_to_csv("GSE6575/CEL")

