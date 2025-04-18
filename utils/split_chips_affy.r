# --- separa_cel.R ------------------------------------------------------------
# Richiede solo affyio (nessuna normalizzazione, oligo/affy non servono qui)
if (!requireNamespace("affyio", quietly = TRUE))
    install.packages("affyio")

library(affyio)

separa_cel <- function(cel_dir) {
  # Trova tutti i .CEL (case‑insensitive)
  cel_files <- list.files(cel_dir, pattern = "\\.CEL$", full.names = TRUE,
                          ignore.case = TRUE)

  if (length(cel_files) == 0) {
    message("Nessun file .CEL trovato in: ", cel_dir)
    return(invisible(NULL))
  }

  # Crea due sottocartelle nella directory sorgente
  affy_dir  <- file.path(cel_dir, "AFFYIO")
  oligo_dir <- file.path(cel_dir, "OLIGO")
  dir.create(affy_dir,  showWarnings = FALSE)
  dir.create(oligo_dir, showWarnings = FALSE)

  # Loop sui file e smista
  for (f in cel_files) {
    chip <- tryCatch(
      affyio::read.celfile.header(f)$cdfName,
      error = function(e) NA_character_
    )

    if (is.na(chip) || !nzchar(chip)) {
      message("Chip type non rilevato per ", basename(f), "; saltato.")
      next
    }

    if (grepl("HuGene|HTA|ST", chip, ignore.case = TRUE)) {
      # compatibile con oligo
      file.copy(f, oligo_dir, overwrite = FALSE)
      message(basename(f), "  -->  OLIGO  (", chip, ")")
    } else {
      # default: compatibile con affy + affyio
      file.copy(f, affy_dir, overwrite = FALSE)
      message(basename(f), "  -->  AFFYIO (", chip, ")")
    }
  }

  message("Smistamento completato.")
  invisible(NULL)
}

# Esempio d’uso
# separa_cel("Datasets/Affymetrix_U133_Plus_v2/GSE6575/CEL")