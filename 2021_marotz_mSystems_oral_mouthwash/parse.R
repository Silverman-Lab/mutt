parse_2021_marotz_mSystems_oral_mouthwash <- function(raw = FALSE) {
  required_pkgs <- c("tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }

  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2021_marotz_mSystems_oral_mouthwash")

  # ----- File paths -----
  counts_zip    <- file.path(local, "2021_marotz_mSystems_oral_mouthwash.RDS.zip")
  metadata_zip  <- file.path(local, "T3_SRS_metadata_ms.txt.zip")

  repro_counts_zips <- c(
    file.path(local, "ERP111447_dada2_counts.rds.zip"),
    file.path(local, "ERP117149_dada2_counts.rds.zip")
  )

  repro_tax_zips <- c(
    file.path(local, "ERP111447_dada2_taxa.rds.zip"),
    file.path(local, "ERP117149_dada2_taxa.rds.zip")
  )

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

  # ----- Original Counts and Taxonomy -----
  orig_rds_file <- unzip(counts_zip, list = TRUE)$Name[1]
  counts_original <- readRDS(unz(counts_zip, orig_rds_file))
  raw_tax <- rownames(counts_original)
  proportions_original <- apply(counts_original, 2, function(col) col / sum(col))

  raw_tax <- data.frame(Taxa = raw_tax)
  tax_original <- raw_tax %>%
    mutate(taxa = str_trim(Taxa)) %>%
    separate(
      Taxa,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = "\\s*;\\s*",
      extra = "drop",
      fill = "right"
    )
  rownames(tax_original) <- rownames(counts_original)
  tax_original = make_taxa_label(tax_original)

  if (!raw) {
    matched_taxa <- tax$Taxa[match(colnames(counts_original), rownames(tax_original))]
    colnames(counts_original) <- matched_taxa
    counts_original <- as.data.frame(t(rowsum(t(counts_original), group = colnames(counts_original))))
  }

  # ----- Metadata and Scale -----
  metadata_txt <- unzip(metadata_zip, list = TRUE)$Name[1]
  metadata <- read.csv(unz(metadata_zip, metadata_txt), sep = "\t", row.names = 1)
  metadata <- metadata[colnames(counts_original), ]
  scale <- metadata$FC_avg_cells_per_ul
  names(scale) <- rownames(metadata)

  # # Process multiple zipped RDS files
  # counts_reprocessed_list <- list()
  # proportions_reprocessed_list <- list()
  # tax_reprocessed_list <- list()

  # for (i in seq_along(repro_counts_zips)) {
  #     # Unzip and read counts
  #     temp_rds <- tempfile(fileext = ".rds")
  #     unzip(repro_counts_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
  #     rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
  #     if (length(rds_files) == 0) stop(paste("No *_counts.rds file found for index", i))
  #     counts <- as.data.frame(readRDS(rds_files[1]))

  #     # Unzip and read taxonomy
  #     temp_tax <- tempfile(fileext = ".rds")
  #     unzip(repro_tax_zips[i], exdir = dirname(temp_tax), overwrite = TRUE)
  #     tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
  #     if (length(tax_files) == 0) stop(paste("No *_taxa.rds file found for index", i))
  #     tax <- as.data.frame(readRDS(tax_files[1]))
  #     tax <- make_taxa_label(tax)

  #     if (!raw) {
  #       matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
  #       colnames(counts) <- matched_taxa
  #       counts <- as.data.frame(t(rowsum(t(counts), group = colnames(counts))))
  #     }

  #     # Compute proportions
  #     proportions <- counts
  #     proportions[] <- lapply(counts, function(col) col / sum(col))

  #     # Label with study name based on zip filename prefix
  #     study_id <- sub("_.*$", "", basename(tools::file_path_sans_ext(repro_counts_zips[i])))
  #     counts$Study <- study_id
  #     proportions$Study <- study_id
  #     tax$Study <- study_id

  #     counts_reprocessed_list[[i]] <- counts
  #     proportions_reprocessed_list[[i]] <- proportions
  #     tax_reprocessed_list[[i]] <- tax
  # }

  # # Combine all
  # counts_reprocessed <- bind_rows(counts_reprocessed_list)
  # proportions_reprocessed <- bind_rows(proportions_reprocessed_list)
  # tax_reprocessed <- bind_rows(tax_reprocessed_list)

  # DELETE LATER WHEN REPROCESS HAS FINALIZED:
  counts_zip = file.path(local, "Marotz_2021_16S.csv.zip")
  counts_reprocessed = read_zipped_table(counts_zip) %>% as.data.frame()
  proportions_reprocessed <- sweep(counts_reprocessed, MARGIN = 1,STATS  = rowSums(counts_reprocessed), FUN = "/")
  ############################################

  if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      proportions_original = fill_na_zero_numeric(proportions_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

  # ----- Return structured list -----
  return(list(
    counts = list(
      original = counts_original,
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = proportions_reprocessed
    ),
    tax = list(
      original = tax_original,
      reprocessed = tax_reprocessed
    ),
    scale = scale,
    metadata = metadata
  ))
}
