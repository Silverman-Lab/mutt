parse_2023_maghini_naturebio_metagenomic <- function() {
<<<<<<< Updated upstream
  required_pkgs <- c("tibble", "tidyverse", "readr")
=======
  required_pkgs <- c("tibble", "tidyverse", "readr", "purrr")
>>>>>>> Stashed changes
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)
  library(readr)
<<<<<<< Updated upstream

  # ----- Local base directory -----
  local <- file.path("2023_maghini_naturebiotechnology_samplemeasurement")

  # ----- File paths -----
  motus_zip      <- file.path(local, "PRJNA940499_motus_merged.tsv.zip")
  metaphlan4_zip <- file.path(local, "PRJNA940499_MetaPhlAn_merged.tsv.zip")

  # ----- Initialize everything as NA -----
  counts_original <- NA
  proportions_original <- NA
  tax_original <- NA
  mOTU3_counts <- NA
  mOTU3_proportions <- NA
  mOTU3_tax <- NA
  MetaPhlAn4_counts <- NA
  MetaPhlAn4_proportions <- NA
  MetaPhlAn4_tax <- NA

  # ----- Read taxonomic counts -----
  tax_files  <- paste0(local, "/bracken_", tax_levels, "_reads.txt")
  tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus")
  raw_counts <- setNames(lapply(tax_files, function(f) {
    read.table(f, header = TRUE, sep = "\t", row.names = 1,
             check.names = FALSE, stringsAsFactors = FALSE)
  }), tax_levels)

  renamed <- imap(raw_counts, function(df, level) {
    orig <- rownames(df)
    new  <- paste0("Taxon_", seq_along(orig))
    list(
      counts = `rownames<-`(df, new),
      tax    = tibble(Taxon = new, OriginalName = orig)
    )
  })
         
  counts_original  <- map(renamed, "counts")
  tax_original     <- map(renamed, "tax")

  # ------- Proportions ----------
  proportions_original <- map(counts_original, function(df) {
    rsums <- rowSums(df)
    prop  <- sweep(df, 1, rsums, FUN = "/")
    prop[is.nan(prop)] <- 0
    return(as_tibble(prop, rownames = "Taxon"))
  })

  # ----- Read and filter qPCR data -----
  plate_list <- list.files(path = base_path, pattern = "qPCR_plate", full.names = TRUE)
  qPCR <- plate_list %>%
    map_dfr(read_csv) %>%
    mutate(ID = gsub("_", "-", SampleName)) %>%
    dplyr::select(Plate, ID, Donor, Condition, PCR_Replicate, Replicate, logCopyNumber, CopyNumber) 

  scale = qPCR %>%  
    dplyr::select(ID, logCopyNumber, CopyNumber)
    group_by(ID) %>%
    summarise(
      mean_CopyNumber     = mean(CopyNumber, na.rm = TRUE),
      sd_CopyNumber       = sd(CopyNumber, na.rm = TRUE),
      mean_logCopyNumber  = mean(logCopyNumber, na.rm = TRUE),
      sd_logCopyNumber    = sd(logCopyNumber, na.rm = TRUE),
      n_reps              = n()
     )

  # ----- Extract and deduplicate metadata -----
  metadata <- qPCR %>%
    dplyr::select(ID, Donor, Condition, Replicate) %>%
    distinct()

  # ----- SRA metadata and merge -----
  metadatafromsra_path <- paste0(base_path, "/SraRunTable (3).csv")
  metadatafromsra <- read.csv(metadatafromsra_path, check.names = FALSE) %>%
    mutate(Library.Name = gsub("(_DNA)", "", Library.Name))

  scale$ID <- gsub("-", "_", scale$ID)
  metadata$ID <- gsub("-", "_", metadata$ID)

  metadata <- merge(
    metadatafromsra, metadata,
    by.x = "Library.Name",
    by.y = "ID",
    all = TRUE
  )

  # ----- mOTU3 Reprocessed -----
  if (file.exists(motus_zip)) {
    motus_files <- unzip(motus_zip, list = TRUE)
    motus_filename <- motus_files$Name[grepl("\\.tsv$", motus_files$Name)][1]
    if (!is.na(motus_filename)) {
        temp_dir <- tempdir()
        unzip(motus_zip, files = motus_filename, exdir = temp_dir, overwrite = TRUE)
        motus_path <- file.path(temp_dir, motus_filename)
        df <- read_tsv(motus_path)
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
        proportions <- apply(df, 2, function(col) col / sum(col))
        tax_df <- data.frame(taxa = rownames(df)) %>%
        mutate(taxa = str_trim(taxa)) %>%
        separate(taxa,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                sep = "\\s*;\\s*", extra = "drop", fill = "right")
        rownames(tax_df) <- rownames(df)

=======
  library(purrr)

  local           <- "2023_maghini_naturebiotechnology_samplemeasurement/"
  motus_zip       <- paste0(local, "PRJNA940499_motus_merged.tsv.zip")
  metaphlan4_zip  <- paste0(local, "PRJNA940499_MetaPhlAn_merged.tsv.zip")

  # ----- Initialize everything as NA -----
  counts_original <- NA
  proportions_original <- NA
  tax_original <- NA
  mOTU3_counts <- NA
  mOTU3_proportions <- NA
  mOTU3_tax <- NA
  MetaPhlAn4_counts <- NA
  MetaPhlAn4_proportions <- NA
  MetaPhlAn4_tax <- NA

  # ----- Read taxonomic counts -----
  tax_files  <- paste0(local, "/bracken_", tax_levels, "_reads.txt")
  tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus")
  raw_counts <- setNames(lapply(tax_files, function(f) {
    read.table(f, header = TRUE, sep = "\t", row.names = 1,
             check.names = FALSE, stringsAsFactors = FALSE)
  }), tax_levels)

  renamed <- imap(raw_counts, function(df, level) {
    orig <- rownames(df)
    new  <- paste0("Taxon_", seq_along(orig))
    list(
      counts = `rownames<-`(df, new),
      tax    = tibble(Taxon = new, OriginalName = orig)
    )
  })
         
  counts_original  <- map(renamed, "counts")
  tax_original     <- map(renamed, "tax")

  # ------- Proportions ----------
  proportions_original <- map(counts_original, function(df) {
    rsums <- rowSums(df)
    prop  <- sweep(df, 1, rsums, FUN = "/")
    prop[is.nan(prop)] <- 0
    return(as_tibble(prop, rownames = "Taxon"))
  })

  # ----- Read and filter qPCR data -----
  plate_list <- list.files(path = base_path, pattern = "qPCR_plate", full.names = TRUE)
  qPCR <- plate_list %>%
    map_dfr(read_csv) %>%
    mutate(ID = gsub("_", "-", SampleName)) %>%
    dplyr::select(Plate, ID, Donor, Condition, PCR_Replicate, Replicate, logCopyNumber, CopyNumber) 

  scale = qPCR %>%  
    dplyr::select(ID, logCopyNumber, CopyNumber)
    group_by(ID) %>%
    summarise(
      mean_CopyNumber     = mean(CopyNumber, na.rm = TRUE),
      sd_CopyNumber       = sd(CopyNumber, na.rm = TRUE),
      mean_logCopyNumber  = mean(logCopyNumber, na.rm = TRUE),
      sd_logCopyNumber    = sd(logCopyNumber, na.rm = TRUE),
      n_reps              = n()
     )

  # ----- Extract and deduplicate metadata -----
  metadata <- qPCR %>%
    dplyr::select(ID, Donor, Condition, Replicate) %>%
    distinct()

  # ----- SRA metadata and merge -----
  metadatafromsra_path <- paste0(base_path, "/SraRunTable (3).csv")
  metadatafromsra <- read.csv(metadatafromsra_path, check.names = FALSE) %>%
    mutate(Library.Name = gsub("(_DNA)", "", Library.Name))

  scale$ID <- gsub("-", "_", scale$ID)
  metadata$ID <- gsub("-", "_", metadata$ID)

  metadata <- merge(
    metadatafromsra, metadata,
    by.x = "Library.Name",
    by.y = "ID",
    all = TRUE
  )

  # ----- mOTU3 Reprocessed -----
  if (file.exists(motus_zip)) {
    motus_files <- unzip(motus_zip, list = TRUE)
    motus_filename <- motus_files$Name[grepl("\\.tsv$", motus_files$Name)][1]
    if (!is.na(motus_filename)) {
        temp_dir <- tempdir()
        unzip(motus_zip, files = motus_filename, exdir = temp_dir, overwrite = TRUE)
        motus_path <- file.path(temp_dir, motus_filename)
        df <- read_tsv(motus_path)
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
        proportions <- apply(df, 2, function(col) col / sum(col))
        tax_df <- data.frame(taxa = rownames(df)) %>%
        mutate(taxa = str_trim(taxa)) %>%
        separate(taxa,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                sep = "\\s*;\\s*", extra = "drop", fill = "right")
        rownames(tax_df) <- rownames(df)

>>>>>>> Stashed changes
        mOTU3_counts <- df
        mOTU3_proportions <- proportions
        mOTU3_tax <- tax_df
    }
  }

  # ----- MetaPhlAn4 Reprocessed -----
  if (file.exists(metaphlan4_zip)) {
    metaphlan4_files <- unzip(metaphlan4_zip, list = TRUE)
    metaphlan4_filename <- metaphlan4_files$Name[grepl("\\.tsv$", metaphlan4_files$Name)][1]
    if (!is.na(metaphlan4_filename)) {
        temp_dir <- tempdir()
        unzip(metaphlan4_zip, files = metaphlan4_filename, exdir = temp_dir, overwrite = TRUE)
        path <- file.path(temp_dir, metaphlan4_filename)
        df <- read_tsv(path)
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
        proportions <- apply(df, 2, function(col) col / sum(col))
        tax_df <- data.frame(taxa = rownames(df)) %>%
        mutate(taxa = str_trim(taxa)) %>%
        separate(taxa,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                sep = "\\s*;\\s*", extra = "drop", fill = "right")
        rownames(tax_df) <- rownames(df)

        MetaPhlAn4_counts <- df
        MetaPhlAn4_proportions <- proportions
        MetaPhlAn4_tax <- tax_df
    }
  }

  # ----- Return -----
  return(list(
    counts = list(
      original = counts_original,
      reprocessed = list(
        mOTU3 = mOTU3_counts,
        MetaPhlAn4 = MetaPhlAn4_counts
      )
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = list(
        mOTU3 = mOTU3_proportions,
        MetaPhlAn4 = MetaPhlAn4_proportions
      )
    ),
    tax = list(
      original = tax_original,
      reprocessed = list(
        mOTU3 = mOTU3_tax,
        MetaPhlAn4 = MetaPhlAn4_tax
      )
    ),
    scale = scale,
    metadata = metadata
  ))
}


