parse_2021_ott_environmentalmicrobiome_riverwater <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  if (!is.logical(raw)) {
    stop("raw must be a logical value")
  }
  if (!is.logical(align)) {
    stop("align must be a logical value")
  }

  library(tibble)
  library(tidyverse)
  library(readr)

  # ----- Local base directory -----
  local <- file.path("2021_ott_environmentalmicrobiome_riverwater")

  # ----- File paths -----
  metadata_zip         <- file.path(local, "metadata.csv.zip")
  sra_zip              <- file.path(local, "SraRunTable (30).csv.zip")
  repro_counts_rds_zip <- file.path(local, "PRJEB42314_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJEB42314_dada2_taxa.rds.zip")

  # ----- Initialize everything as NA -----
  counts_original <- NA
  proportions_original <- NA
  tax_original <- NA
  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA

  # ------ metadata and scale -----
  metadata     <- read_zipped_table(metadata_zip)  # NEED Sample_name from Amelie Ott -- May 5th 2025 follow up with her
  sra          <- read_zipped_table(sra_zip) %>%
                      rename(Accession = run) 

  metadata <- full_join(metadata, sra, by = "Sample_name")
  scale <- metadata %>% dplyr::select("Accession", "Sample_name", "S16_rRNA_ml") %>%
                    mutate(log2_S16_rRNA_ml = ifelse(S16_rRNA_ml > 0, log2(S16_rRNA_ml), NA)) %>%
                    mutate(log10_S16_rRNA_ml = ifelse(S16_rRNA_ml > 0, log10(S16_rRNA_ml), NA))

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
  tax_reprocessed = make_taxa_label(tax_reprocessed)

  # ----- Convert accessions to sample IDs / Sequences to Taxa -----
  if (!raw) {
    aligned_counts <- rename_and_align(counts_reprocessed = counts_reprocessed,
                                      metadata = metadata,
                                      scale = scale,
                                      by_col = "Sample_name",
                                      align = align,
                                      study_name = basename(local))
    counts_reprocessed <- aligned_counts$counts_reprocessed
  }
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

  if (!raw) {
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

  return(list(
    counts      = list(original = NA, reprocessed = counts_reprocessed),
    proportions = list(original = NA, reprocessed = proportions_reprocessed),
    tax         = list(original = NA, reprocessed = tax_reprocessed),
    scale       = scale,
    metadata    = metadata
  ))
}