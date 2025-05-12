parse_2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
      stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
          ". Please install them before running this function.")
  }
  if (!is.logical(align)) {
    stop("align must be a logical value")
  }
  if (!is.logical(raw)) {
    stop("raw must be a logical value")
  }

  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein")

  # ----- File paths -----
  counts_zip         <- NA
  metadata_16s_zip   <- file.path(local, "SraRunTable_16s.csv.zip")
  metadata_18s_zip   <- file.path(local, "SraRunTable_18s.csv.zip")
  scale_zip          <- file.path(local, "16S samples_FCM_Chl.csv.zip")

  repro_counts_zips <- c(
  file.path(local, "PRJNA508514_dada2_counts.rds.zip")#,
  #file.path(local, "PRJNA508517_SILVA_counts.rds.zip")
  )

  repro_tax_zips <- c(
  file.path(local, "PRJNA508514_dada2_taxa.rds.zip")#,
  #file.path(local, "PRJNA508517_SILVA_taxa.rds.zip")
  )

  # ----- Metadata -----
  metadata_16s <- read_zipped_table(metadata_16s_zip, row.names = NULL)
  metadata_16s$source <- "16S"
  metadata_18s <- read_zipped_table(metadata_18s_zip, row.names = NULL)
  metadata_18s$source <- "18S"
  metadata = bind_rows(metadata_16s, metadata_18s) 
  metadata <- metadata %>% rename(Accession = Run, SampleID = sampleid)

  # ----- Scale -----
  scale = read_zipped_table(scale_zip, row.names = NULL)
  mergedwmetadata = scale %>% select(SampleID, Line, Station, `Filtered Seawater Vol [L]`, `SurfChl [mg m-3]`) 
  scale = scale %>% select(SampleID, `FCM [cells/ml]`) %>% 
                        mutate(log2_fc_cells_ml = ifelse(`FCM [cells/ml]` > 0, log2(`FCM [cells/ml]`), NA)) %>%
                        mutate(log10_fc_cells_ml = ifelse(`FCM [cells/ml]` > 0, log10(`FCM [cells/ml]`), NA))  

  metadata = merge(metadata, mergedwmetadata, by = "SampleID")

  # ----- Reprocessed Counts and Taxonomy -----
  repro_labels <- c("16S", "18S")

  counts_reprocessed_list      <- list()
  proportions_reprocessed_list <- list()
  tax_reprocessed_list         <- list()

  if (all(file.exists(repro_counts_zips))) {
    for (i in seq_along(repro_counts_zips)) {
      # ----- Reprocessed counts from RDS ZIP -----
      unzipped = unzip(repro_counts_zips[i], exdir = temp_rds, overwrite = TRUE)
      counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
      if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
      counts_reprocessed <- as.data.frame(readRDS(counts_file))

      # ----- Taxonomy reprocessed -----
      unzipped = unzip(repro_tax_zips[i], exdir = temp_rds, overwrite = TRUE)
      tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
      if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
      tax_reprocessed <- as.data.frame(readRDS(tax_file))
      tax_reprocessed <- make_taxa_label(tax_reprocessed)

      # ----- Match taxa and collapse -----
      if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
        counts_reprocessed <- aligned$reprocessed
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
        original_names <- colnames(counts_reprocessed)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
      }

      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)

      # ----- Proportions -----
      proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')

      # ----- Store -----
      label <- repro_labels[i]
      counts_reprocessed_list[[label]]      <- counts_reprocessed
      proportions_reprocessed_list[[label]] <- proportions_reprocessed
      tax_reprocessed_list[[label]]         <- tax_reprocessed
      cleanup_tempfiles(temp_rds)
    }
  }

  # ----- Return structured list -----
  return(list(
    counts = list(
      original = NA,
      reprocessed = counts_reprocessed_list
    ),
    proportions = list(
      original = NA,
      reprocessed = proportions_reprocessed_list
    ),
    tax = list(
      original = NA,
      reprocessed = tax_reprocessed_list
    ),
    scale = scale,
    metadata = metadata
  ))
}