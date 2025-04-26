parse_2023_Fu_Imeta_wasterwater_pathogens <- function() {
  required_pkgs <- c("tibble", "tidyverse", "R.utils")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2023_Fu_Imeta_wasterwater_pathogens")

    # ----- File paths -----
    scale_meta_prop_zip  <- file.path(local, "data.xlsx.gz")
    mock_scale_zip       <- file.path(local, "sheet2.csv.zip")
    mock_prop_zip        <- file.path(local, "sheet1_pivoted.csv.zip")
    metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
    motus_zip            <- file.path(local, "PRJNA860773_motus_merged.tsv.zip")
    metaphlan4_zip       <- file.path(local, "PRJNA860773_MetaPhlAn_merged.tsv.zip")

  # scale
  unzipped_xlsx <- tempfile(fileext = ".xlsx")
#   R.utils::gunzip(scale_zip, destname = unzipped_xlsx, remove = FALSE, overwrite = TRUE)
#   scale_data <- read_excel(unzipped_xlsx, sheet = 4) 

  # metadata
  csv_name <- unzip(metadata_zip, list = TRUE)$Name[1]
  metadata <- read.csv(unz(metadata_zip, csv_name), check.names = FALSE)
  metadata_1 <- read_excel(unzipped_xlsx, sheet = 2) 

  # proportions
  proportions <- read_excel(unzipped_xlsx, sheet = 5)
  # theres metadata in this

  # tax 
  tax <- tibble(
    Taxon = taxon_names,
    OriginalName = colnames(proportions)["classification"]
  )

  # ----- Initialize everything as NA -----
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

  return(list(
    counts      = list(original = NA, reprocessed = counts_reprocessed),
    proportions = list(original = proportions_original, reprocessed = proportions_reprocessed),
    tax         = list(original = tax, reprocessed = tax_reprocessed),
    scale       = scale,
    metadata    = metadata
  ))
}