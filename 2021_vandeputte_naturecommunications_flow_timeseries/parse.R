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
  #repro_counts_rds_zip<- file.path(local, "EGAS00001005686_dada2_counts.rds.zip")
  #repro_tax_zip       <- file.path(local, "EGAS00001005686_dada2_taxa.rds.zip")
  
  # placeholders
  counts      <- proportions <- tax <- scale <- metadata <- NULL

  # ----- scale, metadata -----
  scale    <- read_zipped_table(scale_zip, sep = "\t", row.names = NULL) %>% 
                mutate(log2_FC_g_mean = ifelse(10^Cell_count_per_gram > 0, log2(10^Cell_count_per_gram), NA)) %>%
                rename(log10_FC_g_mean = Cell_count_per_gram)
  metadata <- read_zipped_table(metadata_zip, row.names = NULL)


  # ----- Read counts (keep first column) -----
  counts_df <- read_zipped_table(counts_zip, sep = "\t", row.names = NULL)
  if (!is.null(counts_df)) {
    id_col          <- counts_df[[1]]            # first column values
    counts_matrix   <- as.matrix(counts_df[ , -1, drop = FALSE])
    rownames(counts_matrix) <- id_col            # give them row names
    counts          <- counts_matrix

    if (!raw) {
      aligned = rename_and_align(counts_original = counts, metadata=metadata, scale=scale, by_col="Sample_ID", align = align, study_name=basename(local))
      counts = aligned$counts_original
    }

    tax      <- read_zipped_table(tax_zip, row.names = NULL)
    tax$Taxa <- tax[[1]]
    tax      <- tax[ , -1] 

    row_sums        <- rowSums(counts, na.rm = TRUE)
    proportions     <- sweep(counts, 1, row_sums, FUN = "/")
    proportions[is.nan(proportions)] <- 0
  } else {
    counts <- proportions <- NULL
  }

  counts_reprocessed = NA
  proportions_reprocessed = NA
  tax_reprocessed = NA

  # # ----- Reprocessed counts from RDS ZIP -----
  # temp_rds <- tempfile(fileext = ".rds")
  # unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

  # rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
  # if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
  # counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

  # # ----- Taxonomy reprocessed -----
  # temp_tax <- tempfile(fileext = ".rds")
  # unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

  # tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
  # if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
  # tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

  
  # # ----- Convert sequences to lowest rank taxonomy found and update key -----
  # tax_reprocessed = make_taxa_label(tax_reprocessed)

  # # ----- Convert accessions to sample IDs / Sequences to Taxa -----
  # if (!raw) { 
  # aligned = rename_and_align(counts_reprocessed = counts, metadata=metadata, scale=scale, by_col="Sample_ID", align = align, study_name=basename(local))
  # counts_reprocessed = aligned$reprocessed
  # }

  # # taxa
  # if (!raw) {
  #     matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
  #     colnames(counts_reprocessed) <- matched_taxa
  #     counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
  # }

  # # proportions reprocessed
  # proportions_reprocessed = counts_reprocessed
  # proportions_reprocessed[-1] <- lapply(
  #     counts_reprocessed[-1],
  #     function(col) col / sum(col)
  # )

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
