parse_2021_barlow_microbiomee_reimaginestudypatientsdpcr <- function(raw = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
        "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
        ". Please install them before running this function."
        )
    }

    library(tidyverse)
    library(readxl)
    library(stringr)
    library(readr)
    
    # ------- local path -----------------
    localPath <- file.path("2021_barlow_microbiomee_reimaginestudypatientsdpcr")

    # ------- file path -----------------
    repro_counts_rds_zip<- file.path(localPath, "PRJNA674353_dada2_counts.rds.zip")
    repro_tax_zip       <- file.path(localPath, "PRJNA674353_dada2_taxa.rds.zip")
    metadata_sra_zip    <- file.path(local, "SraRunTable (33).csv.zip")
    scale_zip           <- file.path(local, "")

    read_zipped_table <- function(zip_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE) {
      if (file.exists(zip_path)) {
        inner_file <- unzip(zip_path, list = TRUE)$Name[1]
        con <- unz(zip_path, inner_file)
        read.table(con, sep = sep, header = header, row.names = row.names, check.names = check.names, stringsAsFactors = FALSE)
      } else {
        warning(paste("File not found:", zip_path))
        return(NA)
      }
    }

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
    
    # ------- original data -------------
    ########## NEEDS SOMEONE TO WORK ON#########
    # ------- counts --------------------
    counts = NA
    # ------- proportions ---------------
    proportions = NA
    # ------- taxa ----------------------
    taxa = NA

    # ------- scale ---------------------
    # ------- metadata ------------------
    sra <- read_zipped_table(metadata_sra_zip)

    # ----- Reprocessed counts from RDS ZIP -----
    temp_rds <- tempfile(fileext = ".rds")
    unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

    rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
    if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

    # ----- Taxonomy reprocessed -----
    temp_tax <- tempfile(fileext = ".rds")
    unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

    tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
    if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

    
    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    make_taxa_label <- function(df) {
        tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
        prefixes  <- c("k", "p", "c", "o", "f", "g")
        if (!all(tax_ranks %in% colnames(df))) {
        stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
        }
        df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
        x[is.na(x) | trimws(x) == ""] <- "unclassified"
        return(x)
        })
        df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
        for (i in length(tax_ranks):1) {
            if (tax_row[i] != "unclassified") {
            return(paste0(prefixes[i], "_", tax_row[i]))
            }
        }
        return("unclassified")  
        })
        return(df)
    }
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    # accessions to sampleIDs is study specific: IF NEED BE

    # taxa
    if (!raw) {
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
    }

    # proportions reprocessed
    proportions_reprocessed = counts_reprocessed
    proportions_reprocessed[-1] <- lapply(
        counts_reprocessed[-1],
        function(col) col / sum(col)
    )
    
    
    return(list(
        scale=scale, 
        metadata=metadata, 
        counts=list(
                    original = counts, 
                    reprocessed = counts_reprocessed
                ),
        proportions=list(
                    original = proportions,
                    reprocessed = proportions_reprocessed
                ),
        tax=list(
                    original = taxa,
                    reprocessed = tax_reprocessed
                )
        )
    )

}