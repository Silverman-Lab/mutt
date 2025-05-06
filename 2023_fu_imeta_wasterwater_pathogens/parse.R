parse_2023_fu_imeta_wasterwater_pathogens <- function(raw = FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "R.utils")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2023_fu_imeta_wasterwater_pathogens")

  # ----- File paths -----
  mock_scale_zip       <- file.path(local, "mock_scale.csv.zip")
  mock_prop_zip        <- file.path(local, "sheet1_pivoted.csv.zip")
  samples_prop_zip     <- file.path(local, "sheet_s5_cleaned.csv.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  motus_zip            <- file.path(local, "PRJNA860773_motus_merged.tsv.zip")
  metaphlan4_zip       <- file.path(local, "PRJNA860773_MetaPhlAn_merged.tsv.zip")

  make_taxa_label <- function(df) {
      tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
      prefixes  <- c("k", "p", "c", "o", "f", "g", "s")
      if (!all(tax_ranks %in% colnames(df))) {
          stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
      }
      df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
          x[is.na(x) | trimws(x) == ""] <- "unclassified"
          x
      })
      df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
          if (tax_row["Genus"] != "unclassified") {
          return(paste0("s_", tax_row["Species"]))
          }
          for (i in (length(tax_ranks)-1):1) { 
          if (tax_row[i] != "unclassified") {
              return(paste0("uc_", prefixes[i], "_", tax_row[i]))
          }
          }
          return("unclassified")
      })
      return(df)
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

  # scale
  mockscale = read_zipped_table(mock_scale_zip, row.names=NULL)

  # metadata
  sra = read_zipped_table(metadata_zip, row.names=NULL)
  # proportions

  proportions = read_zipped_table(samples_prop_zip, row.names=NULL)
  proportions <- proportions %>%
  mutate(`Relative abundance %` = as.numeric(`Relative abundance %`)) %>%
  filter(!is.na(`Relative abundance %`), !is.na(classification), classification != "")  # drop missing classifications

  # Pivot into wide format: rows are bins, columns are classifications
  proportions <- proportions %>%
    select(bins, classification, `Relative abundance %`) %>%
    pivot_wider(
      names_from = classification,
      values_from = `Relative abundance %`,
      values_fill = 0
    ) %>%
    column_to_rownames("bins")

  # tax 
  tax <- tibble(
    taxonomy = colnames(proportions)
  )

  tax$ogtaxonomy <- tax$taxonomy
    tax <- tax %>%
    separate(taxonomy,
            into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
            sep = ";", fill = "right") %>%
    mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))
    tax <- tax %>%
        mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), NA, .)))
    tax[is.na(tax)] <- "unclassified"
    tax = make_taxa_label(tax)
    rownames(tax) <- tax$ogtaxonomy
    if (!raw) {
        matched_taxa <- tax$Taxa[match(colnames(proportions), rownames(tax))]
        colnames(proportions) <- matched_taxa
        proportions <- as.data.frame(t(rowsum(t(proportions), group = colnames(proportions))))
    }

  mockproportions = read_zipped_table(mock_prop_zip, row.names=NULL)

  mocktax <- tibble(
    taxonomy = colnames(mockproportions)
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
    counts      = list(original = list(
                        mocks = NA, 
                        samples = NA
                        ),
                      reprocessed = list(
                        mock = NA,
                        samples = counts_reprocessed
                        )
                      ),
    proportions = list(original = list(
                        mock = mockproportions, 
                        samples = proportions # NEEDS TO BE FIXED
                      ),
                      reprocessed = list(
                        mock = NA,
                        samples = proportions_reprocessed 
                      )
                      ),
    tax         = list(original = list(
                        mock = mocktax, 
                        samples = tax
                        ),
                      reprocessed = list(
                        mock = NA,
                        samples = tax_reprocessed
                        ),
                      ),
    scale       = list(
                    mock = scale,
                    samples = NA #IDK
                      ), 
    metadata    = list(
                    mock = mockmetadata,
                    samples = sra
                      )
  ))
}