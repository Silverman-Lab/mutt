parse_2021_reese_cell_chimpanzee <- function() {
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

    # ---------- local paths --------------------
    local               <- file.path("2021_reese_cell_chimpanzee")

    # ---------- file paths ---------------------
    repro_counts_rds_zip<- file.path(local, "PRJEB39807_dada2_merged_nochim.rds.zip")
    repro_tax_zip       <- file.path(local, "PRJEB39807_dada2_taxonomy_merged.rds.zip")
    metadata_zip        <- file.path(local, "SraRunTable.csv.zip")
    scale_zip           <- file.path(local,"1-s2.0-S0960982220316523-mmc3.xlsx.zip")

    # ---- scale and metadata -------------------
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1] # list file inside zip
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- as.data.frame(read.csv(metadata_path, stringsAsFactors = FALSE))
    rownames(metadata) <- metadata$Sample_name
    metadata <- subset(metadata, select = -Sample_name)
    rownames(metadata) <- gsub("c$", "", rownames(metadata))
    scale_xlsx <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_xlsx, exdir = tempdir(), overwrite = TRUE)
    scale <- as.data.frame(read_excel(scale_path, sheet = "S2 Metadata"))
    rownames(scale) <- scale$`Sample ID`
    scale <- scale[rownames(metadata), , drop = FALSE]
    scale_data <- scale %>%
        select(`Sample ID`, `16S copies per g feces`) %>%
        pivot_wider(names_from = `Sample ID`, values_from = `16S copies per g feces`)
    scale_extra <- scale %>% select(-`Sample ID`, -`16S copies per g feces`)
    scale_extra <- scale_extra[rownames(metadata), , drop = FALSE]
    metadata <- cbind(metadata, scale_extra)

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
    
    # ------ proportions reprocessed -----------
    proportions_reprocessed = counts_reprocessed
    proportions_reprocessed[-1] <- lapply(
        counts_reprocessed[-1],
        function(col) col / sum(col)
    )
    
    # ----- Taxonomy reprocessed ---------------
    temp_tax <- tempfile(fileext = ".rds")
    unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
    tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
    taxonomy_matrix <- readRDS(tax_file)
    rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
    tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
    tax_reprocessed = tax_table

    return(list(
        counts = list(
            original = NA,
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = NA,
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = NA,
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))
}
