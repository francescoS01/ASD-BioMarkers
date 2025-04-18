# funziona bene per vecchi file e non binari

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("affy", quietly = TRUE))
  BiocManager::install("affy")

library(affy)

# Lanciare lo script da dentro ASD-BioMarkers
convert_cel_to_csv <- function(cel_dir_path) {

    root = "Datasets/Affymetrix_U133_Plus_v2/"
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
    cel_files <- list.files(absolute_path, pattern = "\\.CEL$", full.names = TRUE)
    
    if (length(cel_files) == 0) {
        cat("Nessun file .CEL trovato in:", absolute_path, "\n")
        return(NULL)
    }

   
    count = 0
    for (f in cel_files) {
    tryCatch({
        hdr <- read.celfile.header(f)
        cat(basename(f), "->", hdr$Cols, "x", hdr$Rows, "\n")
    }, error = function(e) {
        cat("Errore nella lettura:", basename(f), "\n")
        count <<- count + 1
    })
    }

    good_files <- c()
    for (f in cel_files) {
    tryCatch({
        ReadAffy(filenames = basename(f), celfile.path = absolute_path)
        cat("Leggo il file: ", basename(f),"\n")
        good_files <- c(good_files, f)
    }, error = function(e) {
        cat("File incompatibile o danneggiato:", f, "\n")
    })
    }

    # Carica e normalizza solo i file compatibili
    if (length(good_files) == 0) {
        cat("Nessun file .CEL valido trovato per la normalizzazione.\n")
        return(NULL)
    }
    cat("Normalizzazione di tutti i file .CEL compatibili...\n")
    data <- ReadAffy(filenames = basename(good_files), celfile.path = absolute_path)
    espr <- rma(data)
    matrice <- exprs(espr)

    # Salva un unico CSV con tutti i campioni
    output_path <- file.path(csv_dir_path, "data_all.csv")
    write.csv(matrice, file = output_path)
    if (count > 0){
        cat(count)  
    }       
    cat("Conversione completata.\n")

    # cat(length(cel_files) - length(good_files), "file .CEL esclusi dalla conversione perché incompatibili o danneggiati\n")
}

convert_cel_to_csv("GSE18123/CEL/")