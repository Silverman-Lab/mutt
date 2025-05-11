parse_2024_sternes_frontmicrobiol_IBDppiqPCR <- function(raw = FALSE, align = FALSE) {
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
  
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2024_sternes_frontmicrobiol_IBDppiqPCR")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "PRJNA1120972_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1120972_dada2_taxa.rds.zip")
  metadata_zip         <- file.path(local, "SraRunTable (34).csv.zip")
  scale_zip            <- file.path(local, "scale.csv.zip")

  # --- original counts, proportions, tax ---
  counts <- NA
  proportions <- NA
  tax <- NA
  
  # --- metadata and scale ----
  metadata     <- read_zipped_table(metadata_zip, row.names=NULL) %>% as.data.frame() %>% 
                  rename(Sample_name = sampleid, Accession = Run) %>%
                  mutate(Sample_name = as.character(Sample_name))
  scale_raw    <- read_zipped_table(scale_zip, row.names=NULL) %>% as.data.frame() %>% 
                  rename(Sample_name = ID) %>%  mutate(Sample_name = as.character(Sample_name)) %>% as.data.frame() %>% 
                  mutate(log2_DU_bacterial_load = ifelse(DU_bacterial_load > 0, log2(DU_bacterial_load),NA),
                         log2_TI_bacterial_load = ifelse(TI_bacterial_load > 0, log2(TI_bacterial_load),NA),
                         log2_DNS_bacterial_load = ifelse(DNS_bacterial_load > 0, log2(DNS_bacterial_load),NA)) %>% 
                  mutate(log10_DU_bacterial_load = ifelse(DU_bacterial_load > 0, log10(DU_bacterial_load),NA),
                         log10_TI_bacterial_load = ifelse(TI_bacterial_load > 0, log10(TI_bacterial_load),NA),
                         log10_DNS_bacterial_load = ifelse(DNS_bacterial_load > 0, log10(DNS_bacterial_load),NA))

  meta_part    <- scale_raw %>% select(c("Sample_name", "DX_Groups","AGE","age_binned","Sex", 
  "NC_total_score", "NC_binned", "BMI",	"BMI_binned", "SAGIS_binned", "SAGIS_binned2", "PPI")) 
  scale        = scale_raw %>% select(-c("DX_Groups","AGE","age_binned","Sex", 
  "NC_total_score", "NC_binned", "BMI",	"BMI_binned", "SAGIS_binned", "SAGIS_binned2", "PPI")) 
  metadata     = metadata %>% left_join(meta_part, by = "Sample_name")

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
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name",align = align, study_name=basename(local))
      counts_reprocessed = aligned$reprocessed
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