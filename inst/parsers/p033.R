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
                mutate(log10_FC_g_mean = ifelse(Cell_count_per_gram > 0, log10(Cell_count_per_gram), NA)) %>%
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
  counts_reprocessed2 = NA
  proportions_reprocessed2 = NA
  tax_reprocessed2 = NA

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
    temp_rds <- tempfile("repro")
    dir.create(temp_rds)
    unzipped = unzip(repro_counts_rds_zip, exdir = temp_rds, overwrite = TRUE)
    counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(counts_file))

    # ----- rdp16 -----
    if (!file.exists(file.path(local,"rdp16classified.csv.zip"))) {
      if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
          required_pkgs <- c("dada2", "Biostrings")
          missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
          if (length(missing_pkgs) > 0) {
            stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                ". Please install them before running this function.")
          }
          seqs <- Biostrings::DNAStringSet(colnames(counts_reprocessed))
          rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
          tax_reprocessed2 = make_taxa_label(rdpclassified) 
          write.csv(tax_reprocessed2, file = file.path(local, "rdp16classified.csv"), row.names = TRUE)
        } else {
          stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
      }
       
      } else {
        tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
    }

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
        counts_reprocessed2 = aligned$reprocessed
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        matched_taxa2 <- tax_reprocessed2$Taxa[match(colnames(counts_reprocessed2), rownames(tax_reprocessed2))]
        colnames(counts_reprocessed) <- matched_taxa
        colnames(counts_reprocessed2) <- matched_taxa2
        counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
        counts_reprocessed2 <- collapse_duplicate_columns_exact(counts_reprocessed2)
        original_names <- colnames(counts_reprocessed)
        original_names2 <- colnames(counts_reprocessed2)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
        counts_reprocessed2 <- as.data.frame(lapply(counts_reprocessed2, as.numeric), row.names = rownames(counts_reprocessed2), col.names = original_names2, check.names = FALSE)
        proportions_reprocessed2 <- sweep(counts_reprocessed2, 1, rowSums(counts_reprocessed2), '/')
    }

    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, MARGIN = 1,STATS  = rowSums(counts_reprocessed), FUN = "/")

    cleanup_tempfiles(temp_rds)
  }

  if (!raw) {
      counts = fill_na_zero_numeric(counts)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
      proportions = fill_na_zero_numeric(proportions)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
  }

  # ----- Return structured list -----
  return(list(
      counts = list(
          original = counts,
          reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
      ),
      proportions = list(
          original = proportions,
          reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
      ),
      tax = list(
          original = tax,
          reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
      ),
      scale = scale,
      metadata = metadata
  ))
}
