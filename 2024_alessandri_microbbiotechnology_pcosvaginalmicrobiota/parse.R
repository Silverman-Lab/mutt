parse_2024_alessandri_microbbiotechnology_pcosvaginalmicrobiota <- function() {

  # ----- Required packages -----
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2024_alessandri_microbbiotechnology_pcosvaginalmicrobiota")

  # ----- File paths -----
  scale_zip      <- file.path(local, "MBT2-17-e14540-s002 (1).xlsx.zip")
  metadata_zip   <- file.path(local, "SraRunTable (1).csv.zip")
  motus_zip      <- file.path(local, "PRJNA1075117_motus_merged.tsv.zip")
  metaphlan4_zip <- file.path(local, "PRJNA1075117_MetaPhlAn_merged.tsv.zip")

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


  # ----- Metadata -----
  metadata <- NA
  if (file.exists(metadata_zip)) {
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_con <- unz(metadata_zip, metadata_csv)
    metadata <- read.csv(metadata_con, row.names = "Sample.Name") %>%
      as.data.frame() %>%
      rownames_to_column("Sample")
  }

  # ----- Scale -----
  scale <- NA
  if (file.exists(scale_zip)) {
    scale_files <- unzip(scale_zip, list = TRUE)
    if (nrow(scale_files) > 0) {
      scale_xlsx <- scale_files$Name[1]
      temp_dir <- tempdir()
      unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
      scale_path <- file.path(temp_dir, scale_xlsx)
      sheet7 <- read_xlsx(scale_path, sheet = 7)
      sheet8 <- read_xlsx(scale_path, sheet = 8)
      sheet10 <- read_xlsx(scale_path, sheet = 10)
      sheet11 <- read_xlsx(scale_path, sheet = 11)
      merged <- reduce(list(sheet7, sheet8, sheet10, sheet11), full_join, by = "Sample")

      scale <- merged %>%
        select(Sample, `Bacterial total count`) %>%
        column_to_rownames("Sample")
      additional_metadata <- merged %>% select(-`Bacterial total count`)
      if (!is.null(metadata) && !is.na(metadata)) {
        metadata <- metadata %>%
          full_join(additional_metadata, by = "Sample") %>%
          column_to_rownames("Sample")
      }
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
