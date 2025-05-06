parse_2023_kallastu_research_foodscience_food <- function() {
    required_pkgs <- c("tidyverse", "readxl", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }
    # Load needed libraries
    library(tidyverse)
    library(readxl)
    library(readr)

    # -------- local path -----------
    local                       <- file.path("2023_kallastu_research_foodscience_food")

    # -------- file paths -----------
    metadata_zip                <- file.path(local, "SraRunTable.csv.zip")
    repro_counts_rds_zip        <- file.path(local, "PRJNA861123_dada2_counts.rds.zip")
    repro_tax_zip               <- file.path(local, "PRJNA861123_dada2_taxa.rds.zip")

    # ------- metadata -------------
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- as.data.frame(read.csv(metadata_path, row.names = NULL, stringsAsFactors = FALSE)) %>%
                    rename(Accession = Run, Sample = Sample.Name)

    # ------- scale ----------------
    ## Scale is directly presented as a table in the paper
    scale_data <- data.frame(
        "20St" = c(1.55e9, NA, 1.56e9, 1.16e9),
        `K6 PMA` = c(1.28e8, 1.53e7, 5.96e7, 1.71e8),
        `K6 TOT` = c(1.86e9, NA, 1.26e8, 3.06e8),
        `K5 PMA` = c(3.82e8, 3.64e7, 7.04e7, 3.73e8),
        `K5 TOT` = c(5.90e8, NA, 1.34e8, 5.66e8),
        `K4 PMA` = c(3.97e8, 1.17e8, 2.61e7, 4.12e8),
        `K4 TOT` = c(8.56e8, NA, 1.49e8, 4.60e8),
        `20St 2.5PMA` = c(1.51e9, NA, NA, NA),
        `20St 1PMA`   = c(1.39e9, NA, NA, NA),
        `20St 0.5PMA` = c(1.74e9, NA, NA, NA)
    )
    rownames(scale_data) <- c("Spike_In", "Plating_cfu_per_g", "FC", "qPCR")
    colnames(scale_data) <- c(
        "20St", "K6 PMA", "K6 TOT", "K5 PMA", "K5 TOT",
        "K4 PMA", "K4 TOT", "20St 2.5PMA", "20St 1PMA", "20St 0.5PMA"
    )
    scale = as.data.frame(t(scale_data)) %>% rownames_to_column("Sample")

    # ----- Reprocessed counts from RDS ZIP -----
    temp_rds            <- tempfile(fileext = ".rds")
    unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)
    rds_file            <- list.files(dirname(temp_rds), pattern = "\\.rds$", full.names = TRUE)[1]
    seqtab_nochim       <- readRDS(rds_file)
    rpt_mat             <- t(seqtab_nochim)
    counts_reprocessed  <- as.data.frame(rpt_mat)
    counts_reprocessed$Sequence <- rownames(counts_reprocessed)
    counts_reprocessed = counts_reprocessed[, c("Sequence", setdiff(names(counts_reprocessed), "Sequence"))]
    rownames(counts_reprocessed) <- paste0("Taxon_", seq_len(nrow(counts_reprocessed)))
    
    # proportions reprocessed
    proportions_reprocessed = counts_reprocessed
    proportions_reprocessed[-1] <- lapply(
        counts_reprocessed[-1],
        function(col) col / sum(col)
    )
    
    # ----- Taxonomy reprocessed -----
    temp_tax <- tempfile(fileext = ".rds")
    unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
    tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
    taxonomy_matrix <- readRDS(tax_file)
    rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
    tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
    tax_reprocessed = tax_table
    
    
    return(list(
        scale=scale, 
        metadata=metadata, 
        counts=list(
            original = NA, 
            reprocessed = counts_reprocessed
        ),
        proportions=list(
            original = NA,
            reprocessed = proportions_reprocessed
        ),
        tax=list(
            original = NA,
            reprocessed = tax_reprocessed
        )
        )
    )
}