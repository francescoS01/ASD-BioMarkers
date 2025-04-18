if (!requireNamespace("GEOquery", quietly = TRUE))
    BiocManager::install("GEOquery")

library(GEOquery)
# Increase download timeout and set method for GEOquery
options(timeout = 600)
options(download.file.method.GEOquery = "libcurl")

get_Info_from_GEO <- function(geo, output_dir) {

    # Se il file esiste già, salta l'esecuzione
    output_file <- file.path(output_dir, paste0(geo, "_diagnosis_mapping.csv"))
    if (file.exists(output_file)) {
        cat("File già presente in:", output_file, "\n")
        return(invisible(NULL))
    }

    
    gse <- getGEO(geo, GSEMatrix = TRUE)
    len = length(gse)
    # ci sono casi in cui ci sono più seriex matrix se ci sono più piattaforme, dunque in questo caso aggiunge il nome della piattaforma
    for(i in seq_along(gse)) {
        eset <- gse[[i]]
        if(len > 1){
            # Extract only the platform identifier (text between "-" and "_")
            name_parts <- strsplit(names(gse)[i], "-", fixed = TRUE)[[1]]
            cat(name_parts)
            if (length(name_parts) >= 2) {
                cat("BOIA")
                matrix_name <- strsplit(name_parts[2], "_", fixed = TRUE)[[1]][1]
            } else {
                matrix_name <- names(gse)[i]
            }
        }
        else {
            matrix_name <- NULL
        }


        # estrai la metadata table
        pdata <- pData(eset)

        # Identify all 'characteristics' columns
        char_cols <- grep("^characteristics", colnames(pdata), value = TRUE)

        print(char_cols)
        # Extract diagnosis or fallback to all characteristics
        dx <- sapply(seq_len(nrow(pdata)), function(i) {
            # combine all characteristic strings for sample i
            char_vals <- unlist(strsplit(paste(pdata[i, char_cols], collapse = ";"), ";"))
            char_vals <- trimws(char_vals)
            # find any field starting with 'dx', 'diagnosis', or 'disease'
            dx_entries <- char_vals[grepl("^(dx|diagnosis|disease|Autism trait)", char_vals, ignore.case = TRUE)]
            if (length(dx_entries) == 0) {
                # If no diagnosis-related fields, return all characteristics
                return(paste(char_vals, collapse = "; "))
            } else {
                # remove everything up to the first colon and return the value
                return(sub("^[^:]+:\\s*", "", dx_entries[1]))
            }
        })

        # costruiamo la tabella di mapping
        mapping <- data.frame(
            SampleID = rownames(pdata),
            Title     = as.character(pdata$title),
            Diagnosis_or_Characteristics = dx,
            stringsAsFactors = FALSE
        )

        # Verifica se la directory di output esiste, altrimenti creala
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }

        # Salva il file nella directory specificata
        if (!is.null(matrix_name)){
            output_file <- file.path(output_dir, paste0(geo,"-",matrix_name,"_diagnosis_mapping.csv"))
        }
        else {
            output_file <- file.path(output_dir, paste0(geo,"_diagnosis_mapping.csv"))
        }
        write.csv(mapping, file = output_file, row.names = FALSE)

        cat("File salvato in:", output_file, "\n")
    }

}


if (!interactive() && sys.nframe() == 0){
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Uso: Rscript get_geo_metadata.R <GSE_code> <output_dir>")
}
geo <- args[1]
out_dir <- args[2]
get_Info_from_GEO(geo, out_dir)
}