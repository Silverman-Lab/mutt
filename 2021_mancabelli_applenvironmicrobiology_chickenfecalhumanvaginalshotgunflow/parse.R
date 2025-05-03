parse_2021_mancabelli_applenvironmicrobiology_chickenfecalhumanvaginalshotgunflow <- function(raw = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2021_mancabelli_applenvironmicrobiology_chickenfecalhumanvaginalshotgunflow")

    # ----- File paths -----
    motus_zip            <- file.path(local, "PRJNA641015_motus_merged.tsv.zip")
    metaphlan4_zip       <- file.path(local, "PRJNA641015_MetaPhlAn_merged.tsv.zip")
    metadata_zip         <- file.path(local, "SraRunTable (38).csv.zip")
    scale_zip            <- file.path(local, "")
    proportions_zip      <- file.path(local, "proportions_original.csv.zip")

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

    # ------ original counts ------

    # ------ original proportions and taxa ------
    proportions_original     <- read_zipped_table(proportions_zip)
    tax_original             = as.data.frame(Taxa = proportions_original['Taxonomy'], stringsAsFactors=FALSE)
    proportions_original     = as.data.frame(t(proportions_original))

    # ---- scale and metadata -----
    scale_df     <- read_zipped_table(scale_zip)
    metadata_df  <- read_zipped_table(metadata_zip)

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
    counts = list(
                original = counts_original,
                reprocessed = list(
                        mOTU3 = mOTU3_counts,
                        MetaPhlan4 = MetaPhlAn4_counts
                    )
    ),
    proportions = list(
                original = proportions_original,
                reprocessed = list(
                        mOTU3 = mOTU3_proportions,
                        MetaPhlan4 = MetaPhlAn4_proportions
                    )
    ),
    tax = list(
                original = tax_original,
                reprocessed = list(
                        mOTU3 = mOTU3_tax,
                        MetaPhlan4 = MetaPhlAn4_tax
                    )
    ),
    scale = scale,
    metadata = metadata
    ))
}