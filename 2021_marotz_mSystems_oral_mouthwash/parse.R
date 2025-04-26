parse_2021_marotz_mSystems_oral_mouthwash <- function() {
  required_pkgs <- c("tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }

  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2021_marotz_mSystems_oral_mouthwash")

  # ----- File paths -----
  counts_zip    <- file.path(local, "2021_marotz_mSystems_oral_mouthwash.RDS.zip")
  metadata_zip  <- file.path(local, "T3_SRS_metadata_ms.txt.zip")

  repro_counts_zips <- c(
    file.path(local, "ERP111447_dada2_merged_nochim.rds.zip"),
    file.path(local, "ERP117149_dada2_merged_nochim.rds.zip")
  )

  repro_tax_zips <- c(
    file.path(local, "ERP111447_dada2_taxonomy_merged.rds.zip"),
    file.path(local, "ERP117149_dada2_taxonomy_merged.rds.zip")
  )

  # ----- Original Counts and Taxonomy -----
  orig_rds_file <- unzip(counts_zip, list = TRUE)$Name[1]
  counts_original <- readRDS(unz(counts_zip, orig_rds_file))
  raw_tax <- rownames(counts_original)
  rownames(counts_original) <- paste0("Taxon_", seq_len(nrow(counts_original)))
  proportions_original <- apply(counts_original, 2, function(col) col / sum(col))

  raw_tax <- data.frame(taxa = raw_tax)
  tax_original <- raw_tax %>%
    mutate(taxa = str_trim(taxa)) %>%
    separate(
      taxa,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = "\\s*;\\s*",
      extra = "drop",
      fill = "right"
    )
  rownames(tax_original) <- rownames(counts_original)

  # ----- Metadata and Scale -----
  metadata_txt <- unzip(metadata_zip, list = TRUE)$Name[1]
  metadata <- read.csv(unz(metadata_zip, metadata_txt), sep = "\t", row.names = 1)
  metadata <- metadata[colnames(counts_original), ]
  scale <- metadata$FC_avg_cells_per_ul
  names(scale) <- rownames(metadata)

  # ----- Reprocessed Counts and Taxonomy -----
  counts_reprocessed_list <- list()
  tax_reprocessed_list <- list()

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
    counts_reprocessed_list[[i]] <- counts

    # Taxonomy
    temp_tax <- tempfile(fileext = ".rds")
    unzip(repro_tax_zips[i], exdir = dirname(temp_tax), overwrite = TRUE)
    tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
    taxonomy_matrix <- readRDS(tax_file)
    rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
    tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
    tax_reprocessed_list[[i]] <- tax_table
  }

  # Combine reprocessed
  counts_reprocessed <- bind_rows(counts_reprocessed_list)
  proportions_reprocessed <- counts_reprocessed
  proportions_reprocessed[-1] <- lapply(
    counts_reprocessed[-1],
    function(col) col / sum(col)
  )
  tax_reprocessed <- bind_rows(tax_reprocessed_list)

  # ----- Return structured list -----
  return(list(
    counts = list(
      original = counts_original,
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = proportions_reprocessed
    ),
    tax = list(
      original = tax_original,
      reprocessed = tax_reprocessed
    ),
    scale = scale,
    metadata = metadata
  ))
}
