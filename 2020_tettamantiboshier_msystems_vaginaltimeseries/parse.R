parse_2020_tettamantiboshier_msystems_vaginaltimeseries <- function(raw = FALSE) {
  required_pkgs <- c("tibble", "tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2020_tettamantiboshier_msystems_vaginaltimeseries")

  # ----- File paths -----
  counts_zip           <- file.path(local, "originalcounts.csv.zip")
  scale_zip            <- file.path(local, "scale_connect20250425.csv.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA549339_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA549339_dada2_taxa.rds.zip")

  # ---- helper functions ----
  make_taxa_label <- function(df) {
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    prefixes  <- c("k", "p", "c", "o", "f", "g")
    if (!all(tax_ranks %in% colnames(df))) {
        stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
    }
    df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
        x[is.na(x) | trimws(x) == ""] <- "unclassified"
        x
    })
    df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
        if (tax_row["Genus"] != "unclassified") {
        return(paste0("g_", tax_row["Genus"]))
        }
        for (i in (length(tax_ranks)-1):1) {  # skip Genus
        if (tax_row[i] != "unclassified") {
            return(paste0("uc_", prefixes[i], "_", tax_row[i]))
        }
        }
        return("unclassified")
    })
    return(df)
  }
  fill_na_zero_numeric <- function(x) {
    if (missing(x)) return(NULL)
    if (is.data.frame(x)) {
        x[] <- lapply(x, function(y) if (is.numeric(y)) replace(y, is.na(y), 0) else y)
    } else if (is.matrix(x) && is.numeric(x)) {
        x[is.na(x)] <- 0
    } else if (is.list(x)) {
        x <- lapply(x, fill_na_zero_numeric)
    }
    x
  }
  read_zipped_table <- function(zip_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE) {
      if (file.exists(zip_path)) {
      inner_file <- unzip(zip_path, list = TRUE)$Name[1]
      con <- unz(zip_path, inner_file)
      read.table(con, sep = sep, header = header, row.names = row.names, check.names = check.names, stringsAsFactors = FALSE)
      } else {
      warning(paste("File not found:", zip_path))
      return(NA)
      }
  }

  # ----- Read counts & proportions -----
  if (file.exists(counts_zip)) {
    counts_file <- unzip(counts_zip, list = TRUE)$Name[1]
    counts_con  <- unz(counts_zip, counts_file)
    countsdata  <- read.csv(counts_con) %>% as.data.frame()

    columns_to_drop <- c("Participant", "Hours_In_Study", "Sample_ID")
    counts_original <- countsdata[, !(names(countsdata) %in% columns_to_drop)]
    metadatacols    <- countsdata[, names(countsdata) %in% columns_to_drop]

    tax <- tibble(
      Taxa = colnames(countsdata)[!(colnames(countsdata) %in% columns_to_drop)]
    )

    counts <- bind_cols(metadatacols, counts_original)

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
    # This needs to be cleaned up for the originals because im not sure how to link the sample IDs.
  } else {
    metadata <- NA
  }

 # ----- Reprocessed counts from RDS ZIP -----
  temp_rds <- tempfile(fileext = ".rds")
  unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

  rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
  counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

  # ----- Taxonomy reprocessed -----
  temp_tax <- tempfile(fileext = ".rds")
  unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

  tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
  if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
  tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

  
  # ----- Convert sequences to lowest rank taxonomy found and update key -----
  make_taxa_label <- function(df) {
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    prefixes  <- c("k", "p", "c", "o", "f", "g")
    if (!all(tax_ranks %in% colnames(df))) {
        stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
    }
    df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
        x[is.na(x) | trimws(x) == ""] <- "unclassified"
        x
    })
    df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
        if (tax_row["Genus"] != "unclassified") {
        return(paste0("g_", tax_row["Genus"]))
        }
        for (i in (length(tax_ranks)-1):1) {  # skip Genus
        if (tax_row[i] != "unclassified") {
            return(paste0("uc_", prefixes[i], "_", tax_row[i]))
        }
        }
        return("unclassified")
    })
    return(df)
  }
  tax_reprocessed = make_taxa_label(tax_reprocessed)

  # ----- Convert accessions to sample IDs / Sequences to Taxa -----
  # accessions to sampleIDs is study specific: IF NEED BE

  # taxa
  if (!raw) {
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
  }

  # proportions reprocessed
  proportions_reprocessed = counts_reprocessed
  proportions_reprocessed[-1] <- lapply(
      counts_reprocessed[-1],
      function(col) col / sum(col)
  )

  if (!raw) {
    counts_original = fill_na_zero_numeric(counts_original)
    counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
    proportions_original = fill_na_zero_numeric(proportions_original)
    proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

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
