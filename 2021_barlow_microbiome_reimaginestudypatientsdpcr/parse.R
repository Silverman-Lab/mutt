parse_2021_barlow_microbiomee_reimaginestudypatientsdpcr <- function() {
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
    repro_counts_rds_zip<- file.path(localPath, "PRJNA674353_dada2_merged_nochim.rds.zip")
    repro_tax_zip       <- file.path(localPath, "PRJNA674353_dada2_taxonomy_merged.rds.zip")
    
    
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