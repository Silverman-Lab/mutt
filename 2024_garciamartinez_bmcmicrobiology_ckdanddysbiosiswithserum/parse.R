parse_2024_garciamartinez_bmcmicrobiology_ckdanddysbiosiswithserum <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl")
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
  
  # Load needed libraries 
  library(tidyverse)
  library(readxl)
  
  # ----- Local base directory -----
  local <- file.path("2024_garciamartinez_bmcmicrobiology_ckdanddysbiosiswithserum")

  # ----- File paths -----    
  metadata_zip        <- file.path(local, "SraRunTable (19).csv.zip")
  scale_zip           <- file.path(local, "12866_2024_3590_MOESM1_ESM.xlsx.zip")
  repro_counts_rds_zip<- file.path(local, "PRJEB67373_dada2_counts.rds.zip")
  repro_tax_zip       <- file.path(local, "PRJEB67373_dada2_taxa.rds.zip")

  counts_original = NA
  proportions_original = NA
  tax_original = NA

  # ---- Metadata ----
  metadata = read_zipped_table(metadata_zip, row.names=NULL) %>% rename(Accession = Run, Sample = `Sample Name`) %>% mutate(Sample = as.character(Sample))

  # ---- Scale ----
  scale_files <- unzip(scale_zip, list = TRUE)
  scale_xlsx <- scale_files$Name[1]
  temp_dir <- tempfile("scale")
  dir.create(temp_dir)
  unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
  scale_path <- file.path(temp_dir, scale_xlsx)
  cleanup_tempfiles(temp_dir)
  
  scale_raw <- read_xlsx(scale_path, sheet = 1) %>%
    separate(col = 1, into = c("population", "Sample"), sep = "_", remove = TRUE) %>% mutate(Sample = as.character(Sample))
  metadata <- metadata %>%
    left_join(scale_raw %>% select(Sample, population), by = "Sample") 
  scale <- scale_raw %>%
    select(-population) %>% 
                  mutate(log2_Microbial_load = ifelse(`Microbial_load (no. cells/g)` > 0, log2(`Microbial_load (no. cells/g)`),NA)) %>% 
                  mutate(log10_Microbial_load = ifelse(`Microbial_load (no. cells/g)` > 0, log10(`Microbial_load (no. cells/g)`),NA))

  # ----- Reprocessed counts from RDS ZIP -----
  if (file.exists(repro_counts_rds_zip)) {
    temp_dir <- tempfile("repro")
    dir.create(temp_dir)
    unzipped = unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
    counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(counts_file))

    # ----- Taxonomy reprocessed -----
    unzipped = unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))
    
    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, 
                                  by_col="Sample", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed = collapse_duplicate_columns_exact(counts_reprocessed)
        original_names <- colnames(counts_reprocessed)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
    }

    # proportions reprocessed
    proportions_reprocessed = sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), FUN = "/")
    cleanup_tempfiles(temp_dir)
  }

  if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      proportions_original = fill_na_zero_numeric(proportions_original)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

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
  