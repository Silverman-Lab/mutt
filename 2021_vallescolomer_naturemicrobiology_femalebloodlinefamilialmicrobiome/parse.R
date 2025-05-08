parse_2021_vallescolomer_naturemicrobiology_femalebloodlinefamilialmicrobiome <- function(raw = FALSE, align = FALSE) {
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

    # ---------- local paths --------------------
    local               <- file.path("2021_vallescolomer_naturemicrobiology_femalebloodlinefamilialmicrobiome")

    # ---------- file paths ---------------------
    motus_zip = c(
        file.path(local, "EGAS00001005651_motus_merged.tsv.zip"), # These are not available publicly
        file.path(local, "EGAS00001005649_motus_merged.tsv.zip") # These are not available publicly
    )
    metaphlan4_zip = c(
        file.path(local, "EGAS00001005651_MetaPhlAn_merged.tsv.zip"), # These are not available publicly
        file.path(local, "EGAS00001005649_MetaPhlAn_merged.tsv.zip") # These are not available publicly
    )
    #metadata_zip        <- file.path(local, "SraRunTable.csv.zip") # These are not available publicly
    scale_zip           <- file.path(local,"scalemetadata.csv.zip") # has metadata and scale
    # need to email authors asking for counts with sample IDs that match that datasheet^

    # ------ original counts ------
    counts = NA
    proportions = NA
    tax = NA

    # ---- scale and metadata -----
    scale_df     <- read_zipped_table(scale_zip)

    df <- scale_df %>%
        rename(Sample_name = `Sample ID`)
    scale <- df %>%
        select(Sample_name, `Cell counts (cells/g)`) %>%
        rename(cell_counts = `Cell counts (cells/g)`) %>%
        mutate(log2_cell_counts = ifelse(cell_counts > 0, log2(cell_counts), NA)) %>%
        mutate(log10_cell_counts = ifelse(cell_counts > 0, log10(cell_counts), NA))

    metadata <- df %>%
        select(-`Cell counts (cells/g)`)

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
            if (!raw) {
                aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
                df = aligned$reprocessed
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
                aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
                df = aligned$reprocessed
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
        #counts = fill_na_zero_numeric(counts)
        mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
        #proportions = fill_na_zero_numeric(proportions)
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