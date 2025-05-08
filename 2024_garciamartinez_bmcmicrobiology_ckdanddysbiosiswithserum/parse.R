parse_2024_garciamartinez_bmcmicrobiology_ckdflow <- function(raw = FALSE, align = FALSE) {
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
  metadata = read_zipped_table(metadata_zip, row.names=NULL) %>% rename(Accession = Run)
  
  # ---- Scale ----
  scale_files <- unzip(scale_zip, list = TRUE)
  scale_xlsx <- scale_files$Name[1]
  temp_dir <- tempdir()
  unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
  scale_path <- file.path(temp_dir, scale_xlsx)
  
  scale_raw <- read_xlsx(scale_path, sheet = 1) %>%
    separate(col = 1, into = c("population", "Sample.Name"), sep = "_", remove = TRUE)
  metadata <- metadata %>%
    left_join(scale_raw %>% select(Sample.Name, population), by = "Sample.Name") %>% rename(Sample = Sample.Name)
  scale <- scale_raw %>%
    select(-population) %>% rename(Sample = Sample.Name) %>% 
                  mutate(log2_Microbial_load = ifelse(`Microbial_load (no. cells/g)` > 0, log2(`Microbial_load (no. cells/g)`),NA)) %>% 
                  mutate(log10_Microbial_load = ifelse(`Microbial_load (no. cells/g)` > 0, log10(`Microbial_load (no. cells/g)`),NA))

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
  
  # ----- Convert sequences to lowest rank taxonomy found and update key -----
  tax_reprocessed = make_taxa_label(tax_reprocessed)

  # ----- Convert accessions to sample IDs / Sequences to Taxa -----
  if (!raw) {
  aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, 
                            by_col="Sample", align = align, study_name=basename(local))
  counts_reprocessed = aligned$reprocessed
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
  