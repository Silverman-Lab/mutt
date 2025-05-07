parse_2020_ellegaard_currbiology_honeybeeqPCRshotgunmetagenomics <- function(raw = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2020_ellegaard_currbiology_honeybeeqPCRshotgunmetagenomics")

    # ----- File paths -----
    motus_zip            <- file.path(local, "PRJNA598094_motus_merged.tsv.zip")
    metaphlan4_zip       <- file.path(local, "PRJNA598094_MetaPhlAn_merged.tsv.zip")
    scale_meta_zip       <- file.path(local, "")
    counts_meta_zip      <- file.path(local, "")
    metadata_meta_zip    <- file.path(local, "")


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

    # ----- Initialize everything as NA -----
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
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

    if (!raw) {
        counts = fill_na_zero_numeric(counts)
        mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
        proportions = fill_na_zero_numeric(proportions)
        MetaPhlAn4_counts = fill_na_zero_numeric(MetaPhlAn4_counts)
        mOTU3_proportions = fill_na_zero_numeric(mOTU3_proportions)
        MetaPhlAn4_proportions = fill_na_zero_numeric(MetaPhlAn4_proportions)
    }

    return(list(
        counts = list(
                    original = counts,
                    reprocessed = list(
                            mOTU3 = mOTU3_counts,
                            MetaPhlan4 = MetaPhlAn4_counts
                    )

        ),
        proportions = list(
                    original = proportions,
                    reprocessed = list(
                            mOTU3 = mOTU3_proportions,
                            MetaPhlan4 = MetaPhlAn4_proportions
                        )
        ),
        tax = list(
                    original = tax,
                    reprocessed = list(
                            mOTU3 = mOTU3_tax,
                            MetaPhlan4 = MetaPhlAn4_tax
                        )
        ),
        scale = scale,
        metadata = metadata
    ))
}