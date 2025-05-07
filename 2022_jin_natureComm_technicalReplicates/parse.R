parse_2022_jin_natureComm_technicalReplicates <- function(raw = FALSE) {
  required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      ". Please install them before running this function."
    )
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

  read_xlsx_zip <- function(zipfile,
                            sheet = NULL,
                            skip  = 0,
                            tmp   = tempdir()) {
    if (!file.exists(zipfile)) {
      stop("Zip file not found: ", zipfile)
    }
    zinfo <- utils::unzip(zipfile, list = TRUE)
    print(zinfo$Name)
    xlsx_name <- zinfo$Name[grepl("\\.xlsx$", zinfo$Name)][1]
    if (is.na(xlsx_name)) {
      stop("No .xlsx file found inside ", zipfile)
    }
    extracted <- utils::unzip(zipfile,
                              files     = xlsx_name,
                              exdir     = tmp,
                              overwrite = TRUE,
                              junkpaths = TRUE)
    xlsx_path <- file.path(tmp, basename(xlsx_name))
    if (!file.exists(xlsx_path)) {
      stop("Extraction failed—’", xlsx_path, "’ does not exist")
    }
    dat <- readxl::read_xlsx(xlsx_path,
                            sheet = sheet,
                            skip  = skip)
    return(dat)
  }

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
    matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
    colnames(counts) <- matched_taxa
    counts <- as.data.frame(t(rowsum(t(counts), group = colnames(counts))))
  }

  # ---- scale ----
  scale <- read_xlsx_zip(zipfile = supp8_zip, sheet = "Absolute total abundance", skip = 1)
  
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

  # --- Compute proportions from counts ---
  proportions <- sweep(counts, MARGIN = 1,STATS  = rowSums(counts), FUN = "/")
  
  # Initialize storage lists
  all_counts <- list()
  all_props  <- list()
  all_taxa   <- list()

  for (i in seq_along(repro_counts_zips)) {
    counts_zip <- repro_counts_zips[i]
    tax_zip    <- repro_tax_zips[i]

    # ----- Study prefix -----
    study_prefix <- gsub("_dada2_counts\\.rds\\.zip$", "", basename(counts_zip))

    # ----- Unzip and read counts -----
    temp_rds <- tempfile(fileext = ".rds")
    unzip(counts_zip, exdir = dirname(temp_rds), overwrite = TRUE)
    rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
    if (length(rds_files) == 0) stop("No *_counts.rds file found in: ", counts_zip)
    counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

    # ----- Unzip and read taxonomy -----
    temp_tax <- tempfile(fileext = ".rds")
    unzip(tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
    tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
    if (length(tax_files) == 0) stop("No *_taxa.rds file found in: ", tax_zip)
    tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))
    tax_reprocessed <- make_taxa_label(tax_reprocessed)

    # Taxonomy rownames = ASVs/Features: prefix if needed
    tax_reprocessed$BioProject <- study_prefix
    tax_reprocessed$Sequence <- rownames(tax_reprocessed)

    # Rename counts columns using taxonomy
    if (!raw) {
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), tax_reprocessed$Sequence)]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
    }

    # ----- Proportions -----
    proportions_reprocessed <- counts_reprocessed
    proportions_reprocessed[] <- lapply(
      proportions_reprocessed,
      function(col) col / sum(col)
    )

    # Store results
    all_counts[[i]] <- counts_reprocessed
    all_props[[i]]  <- proportions_reprocessed
    all_taxa[[i]]   <- tax_reprocessed
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


