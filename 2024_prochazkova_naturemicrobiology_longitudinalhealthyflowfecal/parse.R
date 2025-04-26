parse_2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal <- function() {
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "PRJNA1233249_dada2_merged_nochim.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1233249_dada2_taxonomy_merged.rds.zip")
  metadata_zip         <- file.path(local, "SraRunTable.csv.zip")
  scale_zip            <- file.path(local, "41564_2024_1856_MOESM4_ESM_.xlsx.zip")

  counts <- NA
  proportions <- NA
  tax <- NA
  
  zip_list <- unzip(metadata_zip, list = TRUE)
  metadata_csv <- zip_list$Name[1]
  metadata_con <- unz(metadata_zip, metadata_csv)
  metadata <- read.csv(metadata_con, row.names = "person_id") %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  scale_files <- unzip(scale_zip, list = TRUE)
  scale_xlsx <- scale_files$Name[1]
  temp_dir <- tempdir()
  unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
  scale_path <- file.path(temp_dir, scale_xlsx)
  scale_raw <- read_xlsx(scale_path, sheet = 1)
  ml_cols <- grep("^Microbial load", names(scale_raw), value = TRUE)
  ml_mean_sd <- ml_cols[grepl("Mean|SD", ml_cols)]
  meta_cols <- setdiff(names(scale_raw), c("ID", ml_cols))
  meta_part <- scale_raw %>%
    select(ID, all_of(meta_cols))
  scale_part <- scale_raw %>%
    select(ID, all_of(ml_mean_sd))
  metadata <- metadata %>%
    left_join(meta_part, by = "ID") %>%
    column_to_rownames("ID")
  
  scale <- scale_part %>%
    column_to_rownames("ID")

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
