parse_2022_jin_natureComm_technicalReplicates <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      ". Please install them before running this function."
    )
  }
  if (!is.logical(raw) || length(raw) != 1) {
    stop("raw must be a boolean")
  }
  if (!is.logical(align) || length(align) != 1) {
    stop("align must be a boolean")
  }

  library(tidyverse)
  library(readxl)
  library(stringr)
  library(readr)
  
  # -----local path ---------------------------
  local <- file.path("2022_jin_natureComm_technicalReplicates")

  # ---------- file paths ---------------------
  repro_counts_zips <- c(
    file.path(local, "PRJNA639639_dada2_counts.rds.zip"),
    file.path(local, "PRJNA639647_dada2_counts.rds.zip"),
    file.path(local, "PRJNA780331_dada2_counts.rds.zip"),
    file.path(local, "PRJNA780361_dada2_counts.rds.zip")
  )

  repro_tax_zips <- c(
    file.path(local, "PRJNA639639_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA639647_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA780331_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA780361_dada2_taxa.rds.zip")
  )
  supp8_zip           <- file.path(local, "Supplementary Data 8.xlsx.zip")

  sra_zips <- c(
    file.path(local, "SraRunTable (40).csv.zip"),
    file.path(local, "SraRunTable (41).csv.zip"),
    file.path(local, "SraRunTable (42).csv.zip"),
    file.path(local, "SraRunTable (43).csv.zip")
  )

  # Read and combine SRA tables
  sra_tables <- lapply(sra_zips, function(zip) {
    if (file.exists(zip)) {
      read_zipped_table(zip, row.names = NULL) %>% rename(Accession = Run, Sample = `Sample Name`)
    } else {
      NULL
    }
  })
  
  # Combine all tables keeping all columns
  sra <- bind_rows(sra_tables[!sapply(sra_tables, is.null)])

  # Remove ^ characters from Sample column
  sra$Sample <- gsub("\\^", "", sra$Sample)

   # ---- scale ----
  scale <- read_xlsx_zip(zipfile = supp8_zip, sheet = "Absolute total abundance", skip = 1) %>% as.data.frame() %>%
              mutate(log2_copies_mg_mean = ifelse(`Absolute total abudance (copies/mg)` > 0, log2(`Absolute total abudance (copies/mg)`), NA)) %>%
              mutate(log10_copies_mg_mean = ifelse(`Absolute total abudance (copies/mg)` > 0, log10(`Absolute total abudance (copies/mg)`), NA)) %>%
              mutate(log2_copies_mg_sd = ifelse(`Standard deviation (N  = 4 or 5)` > 0, log2(`Standard deviation (N  = 4 or 5)`), NA)) %>%
              mutate(log10_copies_mg_sd = ifelse(`Standard deviation (N  = 4 or 5)` > 0, log10(`Standard deviation (N  = 4 or 5)`), NA))

  # ---- metadata ----
  #extract metadata from sample id
  metadata <- as.data.frame(matrix(NA, nrow = nrow(scale), ncol = 5))
  rownames(metadata) <- scale$Sample
  colnames(metadata) <- c("diet", "subject", "location", "sample_type", "technical_replicate")
  
  #first two letter indicate diet type (VD, VS, CE)
  metadata$diet <- substr(scale$Sample, start = 1, stop = 2)
  
  #third letter indicate subject id within group. Note; subject in different group are different mice (i.e., VDa != VSa)
  metadata$subject <- substr(scale$Sample, start = 3, stop = 3)
  
  #sample collected in two different location (distal, proximal)
  metadata$location <- substr(scale$Sample, start = 4, stop = 7)
  
  #sample measured using either filter-residue (cell-sample) and flow-through (ecDNA-sample)
  sampleid <- rownames(metadata)
  metadata$sample_type <- unlist(lapply(strsplit(sampleid, split = "-"), function(x) x[2]))
  metadata$sample_type[is.na(metadata$sample_type)] <- "cell"
  
  #technical replicate if applicable (not every diet group have technical replicates)
  metadata$technical_replicate <- as.numeric(substr(scale$Sample, start = 8, stop = 8))
  metadata = metadata %>% rownames_to_column(var = "Sample")

  metadata = merge(metadata, sra, by = "Sample")


  # ---- counts and tax and proportions ----
  dat <- read_xlsx_zip(zipfile = supp8_zip,sheet = "Sequencing-determined counts",skip = 1)
  counts <- dat[,2:56]

  tax <- dat[,57:62]
  tax <- tax %>% mutate(across(everything(), 
                  ~ trimws(gsub("\\(.*?\\)", "", .))))

  tax = make_taxa_label(tax)
  rownames(tax) <- dat$`cOTU-ID`
  rownames(counts) <- dat$`cOTU-ID`
  counts = as.data.frame(t(counts))

  if (!raw) {
    aligned = rename_and_align(counts_original = counts, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
    counts = aligned$counts_original
    matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
    colnames(counts) <- matched_taxa
    counts <- collapse_duplicate_columns_exact(counts)
    original_names <- colnames(counts)
    counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
  }

  # --- Compute proportions from counts ---
  proportions <- sweep(counts, MARGIN = 1,STATS  = rowSums(counts), FUN = "/")
  
  # Initialize storage lists
  all_counts <- list()
  all_props  <- list()
  all_taxa   <- list()

  if (all(file.exists(repro_counts_zips)) && all(file.exists(repro_tax_zips))) {
    for (i in seq_along(repro_counts_zips)) {
      counts_zip <- repro_counts_zips[i]
      tax_zip    <- repro_tax_zips[i]

      # ----- Study prefix -----
      study_prefix <- gsub("_dada2_counts\\.rds\\.zip$", "", basename(counts_zip))

      # ----- Unzip and read counts -----
      temp_dir <- tempdir()
      unzip(counts_zip, exdir = temp_dir, overwrite = TRUE)
      rds_files <- list.files(temp_dir, pattern = "_counts\\.rds$", full.names = TRUE)
      if (length(rds_files) == 0) stop("No *_counts.rds file found in: ", counts_zip)
      counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))
      cleanup_tempfiles(temp_dir)

      # ----- Unzip and read taxonomy -----
      temp_dir <- tempdir()
      unzip(tax_zip, exdir = temp_dir, overwrite = TRUE)
      tax_files <- list.files(temp_dir, pattern = "_taxa\\.rds$", full.names = TRUE)
      if (length(tax_files) == 0) stop("No *_taxa.rds file found in: ", tax_zip)
      tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))
      cleanup_tempfiles(temp_dir)

      tax_reprocessed <- make_taxa_label(tax_reprocessed)

      # Taxonomy rownames = ASVs/Features: prefix if needed
      tax_reprocessed$BioProject <- study_prefix
      tax_reprocessed$Sequence <- rownames(tax_reprocessed)

      if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), tax_reprocessed$Sequence)]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
        original_names <- colnames(counts_reprocessed)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
      }

      # ----- Proportions -----
      proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')

      # Store results
      all_counts[[i]] <- counts_reprocessed
      all_props[[i]]  <- proportions_reprocessed
      all_taxa[[i]]   <- tax_reprocessed
    }
  }

  # ----- Merge all dataframes -----
  combined_counts <- do.call(rbind, all_counts)
  combined_props  <- do.call(rbind, all_props)
  combined_taxa   <- bind_rows(all_taxa)

  if (!raw) {
      counts = fill_na_zero_numeric(counts)
      proportions = fill_na_zero_numeric(proportions)
      combined_counts = fill_na_zero_numeric(combined_counts)
      combined_props = fill_na_zero_numeric(combined_props)
  }

  return(list(
    scale=scale, 
    metadata=metadata, 
    counts=list(
      original=counts, 
      reprocessed=combined_counts
    ),
    tax=list(
      original=tax,
      reprocessed=combined_props
    ),
    proportions=list(
      original=proportions,
      reprocessed=combined_taxa
    )
  ))
}


