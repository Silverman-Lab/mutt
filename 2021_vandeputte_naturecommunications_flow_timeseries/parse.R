parse_2021_vandeputte_naturecommunications_flow_timeseries <- function(raw = FALSE) {
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
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

  # Helper: unzip, read the first file, KEEP first column
  read_zip_table <- function(zipfile, sep = "\t", stringsAsFactors = FALSE) {
    if (!file.exists(zipfile)) {
      warning("Archive not found: ", zipfile)
      return(NULL)
    }
    tmpdir <- tempfile("unzip_")
    dir.create(tmpdir)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

    files_in_zip <- unzip(zipfile, exdir = tmpdir)
    if (length(files_in_zip) == 0) {
      warning("Zip is empty: ", zipfile)
      return(NULL)
    }
    read.table(files_in_zip[1],
               header           = TRUE,
               sep              = sep,
               row.names        = NULL,      # <-- keep first column!
               check.names      = FALSE,
               stringsAsFactors = stringsAsFactors)
  }

  # ----- Read counts (keep first column) -----
  counts_df <- read_zip_table(counts_zip, sep = "\t")
  if (!is.null(counts_df)) {
    id_col          <- counts_df[[1]]            # first column values
    counts_matrix   <- as.matrix(counts_df[ , -1, drop = FALSE])
    rownames(counts_matrix) <- id_col            # give them row names
    counts          <- counts_matrix

    row_sums        <- rowSums(counts, na.rm = TRUE)
    proportions     <- sweep(counts, 1, row_sums, FUN = "/")
    proportions[is.nan(proportions)] <- 0
  } else {
    counts <- proportions <- NULL
  }

  # ----- Other tables: tax, scale, metadata -----
  tax      <- read_zip_table(tax_zip,      sep = ",")
  tax$Taxa <- tax[[1]]
  tax      <- tax[ , -1] 
  scale    <- read_zip_table(scale_zip,    sep = "\t")
  metadata <- read_zip_table(metadata_zip, sep = ",")

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
  # make_taxa_label <- function(df) {
  #     tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  #     prefixes  <- c("k", "p", "c", "o", "f", "g")
  #     if (!all(tax_ranks %in% colnames(df))) {
  #     stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
  #     }
  #     df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
  #     x[is.na(x) | trimws(x) == ""] <- "unclassified"
  #     return(x)
  #     })
  #     df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
  #     for (i in length(tax_ranks):1) {
  #         if (tax_row[i] != "unclassified") {
  #         return(paste0(prefixes[i], "_", tax_row[i]))
  #         }
  #     }
  #     return("unclassified")  
  #     })
  #     return(df)
  # }
  # tax_reprocessed = make_taxa_label(tax_reprocessed)

  # # ----- Convert accessions to sample IDs / Sequences to Taxa -----
  # # accessions to sampleIDs is study specific: IF NEED BE

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
