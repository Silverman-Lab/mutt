parse_2024_tunsakul_peerj_aerobicvsanaerobicinhealthyvsobesity <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  if (!is.logical(raw) || length(raw) != 1) {
  stop("`raw` must be a single logical value (TRUE or FALSE)")
  }
  if (!is.logical(align) || length(align) != 1) {
  stop("`align` must be a single logical value (TRUE or FALSE)")
  }

  library(tidyverse)
  
  # ----- Local base directory -----
  local <- file.path("2024_tunsakul_peerj_aerobicvsanaerobicinhealthyvsobesity")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "PRJNA1020208_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1020208_dada2_taxa.rds.zip")
  metadata_zip         <- file.path(local, "SraRunTable (38).csv.zip")
  scale_zip            <- file.path(local, "scale.csv.zip") #MANAN TO FINISH

  # --- original counts, proportions, tax ---
  counts <- NA
  proportions <- NA
  tax <- NA
  
  # --- metadata and scale ----
  metadata     <- read_zipped_table(metadata_zip, row.names=NULL) %>% as.data.frame() %>% 
                  rename(Sample = `Sample Name`, Accession = Run) %>%
                  mutate(Sample = as.character(Sample))
  scale        <- read_zipped_table(scale_zip, row.names=NULL) %>% as.data.frame()
  meta_part    <- scale %>% select(c("Sample", "Condition", "Approximate Mean Load"))
  metadata <- metadata %>%
  left_join(meta_part, by = "Sample") %>%
  mutate(log2_qPCR_copies = ifelse(`Approximate Mean Load` > 0, log2(`Approximate Mean Load`), NA)) %>%
  mutate(log10_qPCR_copies = ifelse(`Approximate Mean Load` > 0, log10(`Approximate Mean Load`), NA))
  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
    temp_dir <- tempfile("repro")
    dir.create(temp_dir)
    unzipped <- unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
    counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(counts_file))

    # ----- Taxonomy reprocessed -----
    unzipped <- unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))

    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, counts_original = counts_original, 
                                              proportions_original = proportions_original, metadata=metadata, scale=scale, by_col="Sample", 
                                              align = align)
      counts_reprocessed = aligned$reprocessed
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
      original_names <- colnames(counts_reprocessed)
      counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
    }
    
    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
    cleanup_tempfiles(temp_dir)
  }

  if (!raw) {
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }
  
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
