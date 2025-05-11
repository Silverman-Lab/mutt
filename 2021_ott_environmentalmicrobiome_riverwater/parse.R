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
  metadata_zip_2       <- file.path(local, "Metadatapaper.csv.zip")
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
  metadata     <- read_zipped_table(metadata_zip, row.names=NULL) %>%
  pivot_longer(
    cols = starts_with("Sequencing_ID(replicate"),
    names_to = "Replicate",
    names_pattern = "Sequencing_ID\\(replicate(\\d)\\)",
    values_to = "Sequencing_ID"
  )
  metadata_2   <- read_zipped_table(metadata_zip_2, row.names=NULL) %>%
  pivot_longer(
    cols = starts_with("Sequencing_ID(replicate"),
    names_to = "Replicate",
    names_pattern = "Sequencing_ID\\(replicate(\\d)\\)",
    values_to = "Sequencing_ID"
  )
  sra          <- read_zipped_table(sra_zip) %>%
                      rename(Accession = run) 

  metadata <- full_join(metadata, sra, by = "Sequencing_ID")
  metadata <- full_join(metadata, metadata_2, by = "Sequencing_ID")
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
                                      by_col = "Sequencing_ID",
                                      align = align,
                                      study_name = basename(local))
      counts_reprocessed <- aligned_counts$counts_reprocessed
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
      original_names <- colnames(counts_reprocessed)
      counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
  }

  # proportions reprocessed
  proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')

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