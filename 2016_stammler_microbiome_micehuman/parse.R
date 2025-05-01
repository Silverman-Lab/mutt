parse_2016_stammler_microbiome_micehuman <- function() {
  required_pkgs <- c("stringr", "tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }

  library(stringr)
  library(tidyverse)

<<<<<<< Updated upstream
  # ----- Local base directory -----
  local <- file.path("2016_stammler_microbiome_micehuman")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "PRJEB11953_dada2_merged_nochim.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJEB11953_dada2_taxonomy_merged.rds.zip")
  scale_16s_zip        <- file.path(local, "Stammler2016_scale.csv.zip")
  counts_16s_zip       <- file.path(local, "Stammler_2016_16S.csv.zip")
  metadata_16s_zip     <- file.path(local, "Stammler_2016_metadata.csv.zip")
=======
  local <- file.path("/2016_stammler_microbiome_micehuman/")
  repro_counts_rds_zip<- paste0(local, "PRJEB11953_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- paste0(local, "PRJEB11953_dada2_taxonomy_merged.rds.zip")
  scale_16s_zip     <- paste0(local, "Stammler2016_scale.csv.zip")
  counts_16s_zip    <- paste0(local, "Stammler_2016_16S.csv.zip")
  metadata_16s_zip  <- paste0(local, "Stammler_2016_metadata.csv.zip") 
>>>>>>> Stashed changes

  read_zipped_csv <- function(zip_path) {
    if (file.exists(zip_path)) {
        csv_file <- unzip(zip_path, list = TRUE)$Name[1]
        read.csv(unz(zip_path, csv_file), row.names = 1, check.names = FALSE)
    } else {
        warning(paste("File not found:", zip_path))
        return(NA)
    }
  }

  # ----- Initialize everything as NA -----
  counts_original_mice <- NA
  proportions_original_mice <- NA
  tax_original_mice <- NA
  counts_original_human <- NA
  proportions_original_human <- NA
  tax_original_human <- NA
  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA

  # ------ original counts ------
  counts_original_mice <- read_zipped_csv(counts_16s_zip)

  if (!is.na(counts_original_mice)[1]) {
<<<<<<< Updated upstream
    original_taxa <- colnames(counts_original_mice)
    taxon_ids <- paste0("Taxon_", seq_len(nrow(counts_original_mice)))
    colnames(counts_original_mice) <- taxon_ids
=======
    original_taxa <- rownames(counts_original_mice)
    taxon_ids <- paste0("Taxon_", seq_len(nrow(counts_original_mice)))
    rownames(counts_original_mice) <- taxon_ids
>>>>>>> Stashed changes

    # Create taxa mapping data frame
    tax_original_mice <- data.frame(
      Taxon = taxon_ids,
      Original_Taxa = original_taxa,
      stringsAsFactors = FALSE
    )

    # ------ proportions from counts ------
<<<<<<< Updated upstream
    proportions_original_mice <- t(counts_original_mice)
=======
    proportions_original_mice <- counts_original_mice
>>>>>>> Stashed changes
    proportions_original_mice[] <- lapply(
      proportions_original_mice,
      function(col) {
        if (is.numeric(col)) {
          total <- sum(col, na.rm = TRUE)
          if (total == 0) return(rep(NA, length(col)))
          return(col / total)
        } else {
          return(col)
        }
      }
    )
  } else {
    proportions_original_mice <- NA
    tax_original_mice <- NA
  }

  # ---- scale and metadata -----
  scale     <- read_zipped_csv(scale_16s_zip)
  metadata  <- read_zipped_csv(metadata_16s_zip)

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

  # ----- Return structured list -----
  return(list(
      counts = list(
          original = counts_original_mice,
          reprocessed = counts_reprocessed
      ),
      proportions = list(
          original = proportions_original_mice,
          reprocessed = proportions_reprocessed
      ),
      tax = list(
          original = tax_original_mice,
          reprocessed = tax_reprocessed
      ),
      scale = scale,
      metadata = metadata
  ))
}



