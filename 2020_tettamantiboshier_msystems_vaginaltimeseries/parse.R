parse_2020_tettamantiboshier_msystems_vaginaltimeseries <- function() {
  required_pkgs <- c("tibble", "tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)

<<<<<<< Updated upstream
  # ----- Local base directory -----
  local <- file.path("2020_tettamantiboshier_msystems_vaginaltimeseries")

  # ----- File paths -----
  counts_zip           <- file.path(local, "originalcounts.csv.zip")
  scale_zip            <- file.path(local, "scale_connect20250425.csv.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA549339_dada2_merged_nochim.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA549339_dada2_taxonomy_merged.rds.zip")

  # ----- Read counts & proportions -----
  if (file.exists(counts_zip)) {
    counts_file <- unzip(counts_zip, list = TRUE)$Name[1]
    counts_con  <- unz(counts_zip, counts_file)
    countsdata  <- read.csv(counts_con) %>% as.data.frame()

    columns_to_drop <- c("Participant", "Hours_In_Study")
    counts_original <- countsdata[, !(names(countsdata) %in% columns_to_drop)]
    metadatacols    <- countsdata[, names(countsdata) %in% columns_to_drop]

    taxon_names <- paste0("Taxon_", seq_len(ncol(counts_original)))
    colnames(counts_original) <- taxon_names

    tax <- tibble(
      Taxon = taxon_names,
      OriginalName = colnames(countsdata)[!(colnames(countsdata) %in% columns_to_drop)]
    )

    counts <- bind_cols(metadatacols, counts_original)

=======
  # ----- File paths -----
  local               <- file.path("/2020_tettamantiboshier_msystems_vaginaltimeseries/")
  counts_zip          <- paste0(local, "originalcounts.csv.zip")
  scale_zip           <- paste0(local, "scale.csv.zip")
  metadata_zip        <- paste0(local, "SraRunTable (30).csv.zip")
  repro_counts_rds_zip<- paste0(local, "PRJNA549339_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- paste0(local, "PRJNA549339_dada2_taxonomy_merged.rds.zip")

  # ----- Read counts & proportions -----
  if (file.exists(counts_zip)) {
    counts_file <- unzip(counts_zip, list = TRUE)$Name[1]
    counts_con  <- unz(counts_zip, counts_file)
    countsdata  <- read.csv(counts_con) %>% as.data.frame()

    columns_to_drop <- c("Participant", "Hours_In_Study")
    counts_original <- countsdata[, !(names(countsdata) %in% columns_to_drop)]
    metadatacols    <- countsdata[, names(countsdata) %in% columns_to_drop]

    taxon_names <- paste0("Taxon_", seq_len(ncol(counts_original)))
    colnames(counts_original) <- taxon_names

    tax <- tibble(
      Taxon = taxon_names,
      OriginalName = colnames(countsdata)[!(colnames(countsdata) %in% columns_to_drop)]
    )

    counts <- bind_cols(metadatacols, counts_original)

>>>>>>> Stashed changes
    row_sums    <- rowSums(counts_original)
    prop_mat    <- sweep(as.matrix(counts_original), 1, row_sums, FUN = "/")
    prop_mat[is.nan(prop_mat)] <- 0
    proportions <- bind_cols(metadatacols, as_tibble(prop_mat))

  } else {
    counts      <- NA
    proportions <- NA
    tax         <- NA
  }

  # ----- Scale data -----
  if (file.exists(scale_zip)) {
    scale_file <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_con  <- unz(scale_zip, scale_file)
    scale  = read.csv(scale_con) %>% as.data.frame()

  } else {
    scale <- NA
  }

  # ----- Metadata -----
  if (file.exists(metadata_zip)) {
    metadata_file <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_con  <- unz(metadata_zip, metadata_file)
    metadata = read.csv(metadata_con) %>% as.data.frame()
    # This needs to be cleaned up because im not sure how to link the sample IDs.
  } else {
    metadata <- NA
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

  # ----- Taxonomy reprocessed -----
  temp_tax <- tempfile(fileext = ".rds")
  unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
  tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
  taxonomy_matrix <- readRDS(tax_file)
  rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
  tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
  tax_reprocessed = tax_table

  # ----- Return combined object -----
  return(list(
    counts      = list(
      original = counts,
      reprocessed = counts_reprocessed
      ),
    proportions = list(
      original = proportions,
      reprocessed = proportions_reprocessed
      ),
    tax         = list(
      original = tax,
      reprocessed = tax_reprocessed
      ),
    scale       = scale,
    metadata    = metadata
  ))
}
