parse_2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein <- function() {
    required_pkgs <- c("tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(tidyverse)

    local               <- "2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein/"
    counts_zip          <- NA
    metadata_16s_zip    <- paste0(local, "SraRunTable_16s.csv.zip")
    metadata_18s_zip    <- paste0(local, "SraRunTable_18s.csv.zip")
    repro_counts_zips   <- c(
        paste0(local, "PRJNA508514_dada2_merged_nochim.rds.zip"),
        paste0(local, "PRJNA508517_SILVA_merged_nochim.rds.zip")
    )
    repro_tax_zips      <- c(
        paste0(local, "PRJNA508514_dada2_taxonomy_merged.rds.zip"),
        paste0(local, "PRJNA508517_SILVA_taxonomy_merged.rds.zip")
    )
    scale_zip           <- paste0 (local, "16S samples_FCM_Chl.csv.zip")

  # ----- Metadata -----

  metadata_txt <- unzip(metadata_16s_zip, list = TRUE)$Name[1] 
  metadata_16s <- read.csv(unz(metadata_16s_zip, metadata_txt), row.names = 1)
  metadata_16s$source <- "16S"
  metadata_txt <- unzip(metadata_18s_zip, list = TRUE)$Name[1]
  metadata_18s <- read.csv(unz(metadata_18s_zip, metadata_txt), row.names = 1)
  metadata_18s$source <- "18S"
  metadata = bind_rows(metadata_16s, metadata_18s)

  # ----- Scale -----
  scale_txt <- unzip(scale_zip, list = TRUE)$Name[1]
  scale = read.csv(unz(scale_zip, scale_txt), row.names = 1)

  # ----- Reprocessed Counts and Taxonomy -----
  repro_labels <- c("16S", "18S")

  counts_reprocessed_list      <- list()
  proportions_reprocessed_list <- list()
  tax_reprocessed_list         <- list()
  
  for (i in seq_along(repro_counts_zips)) {
    # Counts
    temp_rds <- tempfile(fileext = ".rds")
    unzip(repro_counts_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
    rds_file <- list.files(dirname(temp_rds), pattern = "\\.rds$", full.names = TRUE)[1]
    seqtab_nochim <- readRDS(rds_file)
    rpt_mat <- t(seqtab_nochim)
    counts <- as.data.frame(rpt_mat)
    counts$Sequence <- rownames(counts)
    counts <- counts[, c("Sequence", setdiff(names(counts), "Sequence"))]
    rownames(counts) <- paste0("Taxon_", seq_len(nrow(counts)))

    # Proportions
    proportions <- counts
    proportions[-1] <- lapply(proportions[-1], function(col) col / sum(col))

    # Taxonomy
    temp_tax <- tempfile(fileext = ".rds")
    unzip(repro_tax_zips[i], exdir = dirname(temp_tax), overwrite = TRUE)
    tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
    taxonomy_matrix <- readRDS(tax_file)
    rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
    tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")

    # Assign by label (16S or 18S)
    counts_reprocessed_list[[repro_labels[i]]]      <- counts
    proportions_reprocessed_list[[repro_labels[i]]] <- proportions
    tax_reprocessed_list[[repro_labels[i]]]         <- tax_table
  }

  # ----- Return structured list -----
  return(list(
    counts = list(
      original = NA,
      reprocessed = counts_reprocessed_list
    ),
    proportions = list(
      original = NA,
      reprocessed = proportions_reprocessed_list
    ),
    tax = list(
      original = NA,
      reprocessed = tax_reprocessed_list
    ),
    scale = scale,
    metadata = metadata
  ))
}