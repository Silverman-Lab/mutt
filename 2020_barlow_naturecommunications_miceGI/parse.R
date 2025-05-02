parse_2020_barlow_naturecommunications_miceGI <- function() {
    required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }

    # Load libraries
    library(tidyverse)
    library(readxl)
    library(stringr)
    library(readr)

    # ----- Local base directory -----
    local <- file.path("2020_barlow_naturecommunications_miceGI")

    # ----- File paths -----
    metadata_zip         <- file.path(local, "metadata.xlsx.zip")
    sra_metadata_zip     <- file.path(local, "SraRunTable (32).csv.zip")
    counts_zip           <- NA  # No original counts provided
    proportions_zip      <- file.path(local, "Relative_Abundance_Table.csv.zip")
    scale_zip            <- file.path(local, "loaddata.csv.zip")
    repro_counts_rds_zip <- file.path(local, "PRJNA575097_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA575097_dada2_taxa.rds.zip")

    # --- proportions ---
    prop_csv <- unzip(proportions_zip, list = TRUE)$Name[1]
    prop_path <- unzip(proportions_zip, files = prop_csv, exdir = tempdir(), overwrite = TRUE)
    proportions <- read.csv(prop_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

    # --- Assign Taxon_1 ... Taxon_N as rownames and store original names ---
    original_taxa <- rownames(proportions)
    taxon_ids <- paste0("Taxon_", seq_len(nrow(proportions)))
    rownames(proportions) <- taxon_ids

    # --- Create taxa dataframe (lookup table) ---
    tax <- data.frame(
    Taxon = taxon_ids,
    Original_Taxa = original_taxa,
    stringsAsFactors = FALSE
    )

    # ----- metadata ------
    metadata_xlsx <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_path <- unzip(metadata_zip, files = metadata_xlsx, exdir = tempdir(), overwrite = TRUE)
    metadata <- read_excel(metadata_path, sheet="Sheet1")

    sra_metadata_csv <- unzip(sra_metadata_zip, list = TRUE)$Name[1]
    sra_metadata_path <- unzip(sra_metadata_zip, files = sra_metadata_csv, exdir = tempdir(), overwrite = TRUE)
    sra_metadata <- read.csv(sra_metadata_path, row.names = 1, stringsAsFactors = FALSE)
    rownames(sra_metadata) <- sra_metadata$`Sample Name`
    sra_metadata <- subset(sra_metadata, select = -Sample.Name)

    # Sample names are all screwed up across SRA, loaddata, and metadata -- cant merge them. 
    

    # --- Scale (qPCR) ---
    scale_csv <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_csv, exdir = tempdir(), overwrite = TRUE)
    scale <- read.csv(scale_path, row.names = 1, stringsAsFactors = FALSE)

    metadata <- metadata[rownames(metadata) %in% colnames(counts), ]
    scale <- scale[rownames(scale) %in% colnames(counts), ]
    sra_metadata <- sra_metadata[rownames(sra_metadata) %in% colnames(counts), ]

    # Merge sra_metadata and metadata
    metadata <- cbind(metadata, sra_metadata[rownames(metadata), , drop = FALSE])

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

    # Return structured list
    return(list(
        counts = list(
            original = counts,
            reprocessed = counts_reprocessed,
        )
        tax = list(
            original = tax,
            reprocessed = tax_reprocessed,
        )
        proportions = list(
            original = proportions,
            reprocessed = proportions_reprocessed,
        )
        metadata = metadata,
        scale = scale,
    ))
}
