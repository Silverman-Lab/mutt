parse_2022_suriano_aps_micefecal <- function() {
  required_pkgs <- c("tidyverse", "readxl", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tidyverse)
  library(readxl)
  library(readr)

  # ----- local path ------------------------
  local                   <- file.path("2022_suriano_aps_micefecal")

  # ----- file paths -------------------------
  metadata_zip            <- file.path(local,"SraRunTable.csv.zip")
  scale_zip               <- file.path(local,"Supplementary Table 6.xlsx.zip")
  repro_counts_rds_zip    <- file.path(local,"PRJEB53668_dada2_merged_nochim.rds.zip")
  repro_tax_zip           <- file.path(local,"PRJEB53668_dada2_taxonomy_merged.rds.zip")

  # ----- metadata ---------------------------
  metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]  # list file inside zip
  metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
  metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

  metadata$Sample_name <- gsub("^(.*?_D)(\\d{2})\\.(\\d+)$", "D\\2-\\3", metadata$Sample_name)
  rownames(metadata) <- metadata$Sample_name
  metadata <- subset(metadata, select = -Sample_name)

  # ----- scale ------------------------------
  scale_xlsx <- unzip(scale_zip, list = TRUE)$Name[1]
  scale_path <- unzip(scale_zip, files = scale_xlsx, exdir = tempdir(), overwrite = TRUE)
  
  scale <- read_excel(scale_path, sheet = "Microbial loads")
  scale <- scale %>%
    select(Sample, `Cells.g.of.fecal.sample`) %>%
    pivot_wider(names_from = Sample, values_from = `Cells.g.of.fecal.sample`)
  
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