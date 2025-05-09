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
  metadata_txt <- unzip(metadata_16s_zip, list = TRUE)$Name[1] 
  metadata_16s <- read.csv(unz(metadata_16s_zip, metadata_txt), row.names = 1)
  metadata_16s$source <- "16S"
  metadata_txt <- unzip(metadata_18s_zip, list = TRUE)$Name[1]
  metadata_18s <- read.csv(unz(metadata_18s_zip, metadata_txt), row.names = 1)
  metadata_18s$source <- "18S"
  metadata = bind_rows(metadata_16s, metadata_18s) %>% rename(Accession = Run)

  # ----- Scale -----
  scale_txt <- unzip(scale_zip, list = TRUE)$Name[1]
  scale = read.csv(unz(scale_zip, scale_txt), row.names = 1)
  mergedwmetadata = scale %>% select(sampleID, Line, Station, `Filtered Seawater Vol [L]`, `SurfChl [mg m-3]`) 
  scale = scale %>% select(sampleID, `FCM [cells/ml]`) %>% 
                        mutate(log2_fc_cells_ml = ifelse(`FCM [cells/ml]` > 0, log2(`FCM [cells/ml]`), NA)) %>%
                        mutate(log10_fc_cells_ml = ifelse(`FCM [cells/ml]` > 0, log10(`FCM [cells/ml]`), NA))  

  metadata = merge(metadata, mergedwmetadata, by.x = "ID", by.y = "sampleID")

  # ----- Reprocessed Counts and Taxonomy -----
  repro_labels <- c("16S", "18S")

  counts_reprocessed_list      <- list()
  proportions_reprocessed_list <- list()
  tax_reprocessed_list         <- list()

  if (all(file.exists(repro_counts_zips))) {
    for (i in seq_along(repro_counts_zips)) {
      # ----- Reprocessed counts from RDS ZIP -----
      temp_rds <- tempfile(fileext = ".rds")
      unzip(repro_counts_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)

      rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
      if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
      counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

      # ----- Taxonomy reprocessed -----
      temp_tax <- tempfile(fileext = ".rds")
      unzip(repro_tax_zips[i], exdir = dirname(temp_tax), overwrite = TRUE)

      tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
      if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
      tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))
      tax_reprocessed <- make_taxa_label(tax_reprocessed)

      # ----- Match taxa and collapse -----
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))

      if (!raw) {
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        align = rename_and_align(counts_reprocessed = counts_reprocessed, metadata, scale, by_col = "sampleID", align = align, study_name = basename(local))
        counts_reprocessed <- align$reprocessed
      }


      # ----- Proportions -----
      proportions_reprocessed <- counts_reprocessed
      proportions_reprocessed[] <- lapply(proportions_reprocessed, function(col) col / sum(col))

      # ----- Store -----
      label <- repro_labels[i]
      counts_reprocessed_list[[label]]      <- counts_reprocessed
      proportions_reprocessed_list[[label]] <- proportions_reprocessed
      tax_reprocessed_list[[label]]         <- tax_reprocessed
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