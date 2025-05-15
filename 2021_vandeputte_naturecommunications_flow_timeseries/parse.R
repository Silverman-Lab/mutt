parse_2021_vandeputte_naturecommunications_flow_timeseries <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl")
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

  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2021_vandeputte_naturecommunications_flow_timeseries")

  # ----- File paths -----
  counts_zip          <- file.path(local, "Vandeputte_2021_16S.tsv.zip")
  scale_zip           <- file.path(local, "Vandeputte_2021_load.tsv.zip")
  tax_zip             <- file.path(local, "tax.csv.zip")
  metadata_zip        <- file.path(local, "metadata.csv.zip")
  #sra_metadata_zip    <- file.path(local, "SraRunTable (38).csv.zip")
  repro_counts_rds_zip<- file.path(local, "EGAS00001005686_dada2_counts.rds.zip")
  repro_tax_zip       <- file.path(local, "EGAS00001005686_dada2_taxa.rds.zip")
  
  # placeholders
  counts      <- proportions <- tax <- scale <- metadata <- NULL
  # ----- scale, metadata -----
  scale    <- read_zipped_table(scale_zip, sep = "\t", row.names = NULL) %>% 
                mutate(log2_FC_g_mean = ifelse(Cell_count_per_gram > 0, log2(Cell_count_per_gram), NA)) %>%
                rename(log10_FC_g_mean = ifelse(Cell_count_per_gram > 0, log10(Cell_count_per_gram), NA)) %>%
                select(-Subject_ID)
  metadata <- read_zipped_table(metadata_zip, row.names = NULL)


  # ----- Read counts -----
  counts <- read_zipped_table(counts_zip, sep = "\t", row.names = 1)
  if (!is.null(counts)) {
    counts = fill_na_zero_numeric(counts)

    tax      <- read_zipped_table(tax_zip, row.names = NULL)
    tax$Taxa <- tax[[1]]
    tax      <- tax[ , -1] 

    if (!raw) {
      aligned = rename_and_align(counts_original = counts, metadata=metadata, scale=scale, by_col="Sample_ID", align = align, study_name=basename(local))
      counts = aligned$counts_original
      original_names <- colnames(counts)
      counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
    }

    row_sums        <- rowSums(counts, na.rm = TRUE)
    proportions     <- sweep(counts, 1, row_sums, FUN = "/")
  } else {
    counts <- proportions <- NULL
  }

  counts_reprocessed = NA
  proportions_reprocessed = NA
  tax_reprocessed = NA

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
    temp_rds <- tempfile("repro")
    dir.create(temp_rds)
    unzipped = unzip(repro_counts_rds_zip, exdir = temp_rds, overwrite = TRUE)
    counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(counts_file))

    # ----- Taxonomy reprocessed -----
    unzipped = unzip(repro_tax_zip, exdir = temp_rds, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))
    
    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) { 
        aligned = rename_and_align(counts_reprocessed = counts, metadata=metadata, scale=scale, by_col="Sample_ID", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
        original_names <- colnames(counts_reprocessed)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
    }

    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, MARGIN = 1,STATS  = rowSums(counts_reprocessed), FUN = "/")

    cleanup_tempfiles(temp_rds)
  }

  if (!raw) {
      counts = fill_na_zero_numeric(counts)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions = fill_na_zero_numeric(proportions)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

  # ----- Return structured list -----
  return(list(
      counts = list(
          original = counts,
          reprocessed = counts_reprocessed
      ),
      proportions = list(
          original = proportions,
          reprocessed = proportions_reprocessed
      ),
      tax = list(
          original = tax,
          reprocessed = tax_reprocessed
      ),
      scale = scale,
      metadata = metadata
  ))
}
