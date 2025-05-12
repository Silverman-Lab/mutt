parse_2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal <- function(raw = FALSE, align = FALSE) {
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
  local <- file.path("2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "PRJNA1233249_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1233249_dada2_taxa.rds.zip")
  metadata_zip         <- file.path(local, "SraRunTable.csv.zip")
  scale_zip            <- file.path(local, "scale.csv.zip")

  # --- original counts, proportions, tax ---
  counts <- NA
  proportions <- NA
  tax <- NA
  
  # --- metadata and scale ----
  metadata     <- read_zipped_table(metadata_zip, row.names=NULL) %>% as.data.frame() %>% 
                  rename(ID = person_id, Accession = Run) %>%
                  mutate(ID = as.character(ID))
  scale_raw    <- read_zipped_table(scale_zip, row.names=NULL) %>% as.data.frame()
  ml_cols      <- grep("^Microbial_load", names(scale_raw), value = TRUE)
  ml_mean_sd   <- ml_cols[grepl("Mean|SD", ml_cols)]
  meta_cols    <- setdiff(names(scale_raw), c("ID", ml_cols))
  meta_part    <- scale_raw %>% select(ID, all_of(meta_cols)) %>% mutate(ID = as.character(ID))
  scale        = scale_raw %>% select(ID, all_of(ml_mean_sd)) %>% mutate(ID = as.character(ID)) %>% 
                  mutate(log2_Microbial_load_Mean = ifelse(Microbial_load_Mean > 0, log2(Microbial_load_Mean),NA),
                         log2_Microbial_load_SD = ifelse(Microbial_load_SD > 0, log2(Microbial_load_SD),NA)) %>% 
                  mutate(log10_Microbial_load_Mean = ifelse(Microbial_load_Mean > 0, log10(Microbial_load_Mean),NA),
                         log10_Microbial_load_SD = ifelse(Microbial_load_SD > 0, log10(Microbial_load_SD),NA))
                        		
  metadata     = metadata %>% left_join(meta_part, by = "ID")

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
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
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="ID",align = align, study_name=basename(local))
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
