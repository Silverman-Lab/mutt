parse_2017_vandeputte_nature_flow <- function() {
  required_pkgs <- c("tibble", "tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2017_vandeputte_nature_flow")

  # ----- File paths -----
  metadata_zip         <- file.path(local, "Vandeputte_2017_metadata.csv.zip")
  metadata_two_zip     <- file.path(local, "cellcountstotal.csv.zip")
  orig_counts_zip      <- file.path(local, "OTU_nochim.zip")
  orig_tax_rdp_zip     <- file.path(local, "otu_taxonomy_rdp.csv.zip")
  orig_tax_silva_zip   <- file.path(local, "otu_taxonomy_silva.csv.zip")
  orig_prop_zip        <- NA
  repro_counts_rds_zip <- file.path(local, "PRJEB21504_dada2_merged_nochim.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJEB21504_dada2_taxonomy_merged.rds.zip")


  # ----- Metadata and Scale -----
  # Read scale
  count_csv    <- unzip(metadata_two_zip, list = TRUE)$Name[1]
  count_con    <- unz(metadata_two_zip, count_csv)
  metadata_two <- read.csv(count_con) %>% as.data.frame()

  # Read metadata
  meta_csv     <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_zip, meta_csv)
  metadata_df  <- read.csv(meta_con) %>% as.data.frame()

  # Join and select
  df <- left_join(
    metadata_df,
    metadata_two %>% select(Sample, SampleID) %>% rename(Accession = SampleID),
    by = c("Individual" = "Sample")
  )
  metadata = df %>% select(
    Individual, Cohort, Day,
    `Health status`, accession = Accession, Enterotype
  )
  scale = df %>% select(
    Individual, Day, accession = Accession,
    `Average cell count (per gram of fresh feces)`,
    `STDEV cell count (per gram of fresh feces)`,
    `Average cell count (per gram of frozen feces)`,
    `STDEV cell count (per gram of frozen feces)`
  )

  # ----- Original counts from CSV.zip -----
  if (file.exists(orig_counts_zip)) {
    orig_csv <- unzip(orig_counts_zip, list = TRUE)$Name[1]
    orig_con <- unz(orig_counts_zip, orig_csv)
    orig_mat <- read.csv(orig_con, row.names = 1, check.names = FALSE)
    counts_original <- as.data.frame(orig_mat)
    counts_original$Sequence <- rownames(counts_original)
    counts_original <- counts_original[, c("Sequence", setdiff(names(counts_original), "Sequence"))]
    rownames(counts_original) <- paste0("Taxon_", seq_len(nrow(counts_original)))
  } else {
    counts_original <- NA
  }

  if (file.exists(orig_counts_zip)) {
    proportions_original = counts_original
    proportions_original[-1] <- lapply(
      counts_original[-1],
      function(col) col / sum(col)
    )
  } else if (file.exists(orig_prop_zip)) {
    prop_csv <- unzip(orig_prop_zip, list = TRUE)$Name[1]
    prop_con <- unz(orig_prop_zip, prop_csv)
    proportions_original = read.csv(prop_con, row.names = 1, check.names = FALSE) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sequence") %>%
      dplyr::select(Sequence, everything())
  } else {
  proportions_original <- NA
  }

  # ----- Reprocessed counts from RDS ZIP -----
  temp_rds            <- tempfile(fileext = ".rds")
  unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)
  rds_file            <- list.files(dirname(temp_rds), pattern = "\\.rds$", full.names = TRUE)[1]
  seqtab_nochim       <- readRDS(rds_file)
  rpt_mat             <- t(seqtab_nochim)
  counts_reprocessed  <- as.data.frame(rpt_mat)
  counts_reprocessed$Sequence <- rownames(counts_reprocessed)
  counts_reprocessed = counts_reprocessed[, c("Sequence", setdiff(names(counts_reprocessed), "Sequence"))]
  rownames(counts_reprocessed) <- paste0("Taxon_", seq_len(nrow(counts_reprocessed)))

  # proportions reprocessed
  proportions_reprocessed = counts_reprocessed
  proportions_reprocessed[-1] <- lapply(
    counts_reprocessed[-1],
    function(col) col / sum(col)
  )

  # --- Original taxonomies ---
  read_taxonomy_zip <- function(zip_path) {
    if (!file.exists(zip_path)) return(NA)
    zip_contents <- unzip(zip_path, list = TRUE)
    csv_name <- zip_contents$Name[1]
    con <- unz(zip_path, csv_name)
    tax_df <- read.csv(con, row.names = 1, check.names = FALSE)
    tax_df$Sequence <- rownames(tax_df)
    tax_df <- tax_df[, c("Sequence", setdiff(names(tax_df), "Sequence"))]
    rownames(tax_df) <- paste0("Taxon_", seq_len(nrow(tax_df)))
    as_tibble(tax_df, rownames = "Taxon")
  }

  tax_original_rdp   <- read_taxonomy_zip(orig_tax_rdp_zip)
  tax_original_silva <- read_taxonomy_zip(orig_tax_silva_zip)

  # ----- Taxonomy reprocessed -----
  temp_tax <- tempfile(fileext = ".rds")
  unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
  tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
  taxonomy_matrix <- readRDS(tax_file)
  rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
  tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
  tax_reprocessed <- tax_table

  # ----- Return all -----
  return(list(
    counts      = list(
      original    = counts_original,
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original    = proportions_original,
      reprocessed = proportions_reprocessed
    ),
    tax         = list(
      original_rdp    = tax_original_rdp,
      original_silva  = tax_original_silva,
      reprocessed = tax_reprocessed
    ),
    scale       = scale,
    metadata    = metadata
  ))
}
