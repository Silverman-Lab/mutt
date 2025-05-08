parse_2013_reyes_procnationalacademyscience_gnotobioticmousehumangutflow <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }
    if (!is.logical(raw)) {
        stop("raw must be a logical value")
    }
    if (!is.logical(align)) {
        stop("align must be a logical value")
    }
    # Load needed libraries
    library(tidyverse)
    library(readxl)
    library(readr)

    # -------- local path -----------
    local                       <- file.path("2013_reyes_procnationalacademyscience_gnotobioticmousehumangutflow")


    # -------- file paths -----------
    scale_zip      <- file.path(local, "sd01.xlsx.zip")
    metadata_zip   <- file.path(local, "SraRunTable (33).csv.zip")
    motus_zip      <- file.path(local, "PRJNA1075117_motus_merged.tsv.zip")
    metaphlan4_zip <- file.path(local, "PRJNA1075117_MetaPhlAn_merged.tsv.zip")

     # ----- Metadata -----
    metadata <- NA
    if (file.exists(metadata_zip)) {
        metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
        metadata_con <- unz(metadata_zip, metadata_csv)
        metadata <- read.csv(metadata_con, row.names = "Sample_name") %>%
        as.data.frame() %>%
        rownames_to_column("Sample")
    }

    # ----- Scale ----- # Needs to be fixed. Supp table 3 and 4 so sheets 3 and 4
    scale <- NA
    if (file.exists(scale_zip)) {
        scale_csv <- unzip(scale_zip, list = TRUE)$Name[1]
        scale_con <- unz(scale_zip, scale_csv)
        scale <- read.csv(scale_con, row.names = "Sample ID") %>%
        as.data.frame() %>%
        rownames_to_column("Sample")
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
            if (!raw) {
                align <- rename_and_align(counts_reprocessed = df, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
                df = align$reprocessed
            }
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
            if (!raw) {
                align <- rename_and_align(counts_reprocessed = df, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
                df = align$reprocessed
            }
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
        counts_original = fill_na_zero_numeric(counts_original)
        mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
        proportions_original = fill_na_zero_numeric(proportions_original)
        MetaPhlAn4_counts = fill_na_zero_numeric(MetaPhlAn4_counts)
        mOTU3_proportions = fill_na_zero_numeric(mOTU3_proportions)
        MetaPhlAn4_proportions = fill_na_zero_numeric(MetaPhlAn4_proportions)
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