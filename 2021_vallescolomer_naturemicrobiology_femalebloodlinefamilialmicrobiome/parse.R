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
    counts = NA
    proportions = NA
    tax = NA
    mOTU3_counts <- NA
    mOTU3_proportions <- NA
    mOTU3_tax <- NA
    MetaPhlAn4_counts <- NA
    MetaPhlAn4_proportions <- NA
    MetaPhlAn4_tax <- NA

    # ----- mOTU3 Reprocessed -----
    if (file.exists(motus_zip)) {
    # 1. create a private scratch folder
    temp_dir <- tempfile("motus_")
    dir.create(temp_dir)

    # 2. find the .tsv inside the ZIP
    motus_files    <- unzip(motus_zip, list = TRUE)
    motus_filename <- motus_files$Name[
                        grepl("\\.tsv$", motus_files$Name, ignore.case = TRUE)
                    ][1]

    if (!is.na(motus_filename)) {
        # 3. extract just that file and grab its full path
        unzipped    <- unzip(
                        motus_zip,
                        files     = motus_filename,
                        exdir     = temp_dir,
                        overwrite = TRUE
                    )
        motus_path  <- unzipped[1]

        # 4. read counts + set rownames
        df <- readr::read_tsv(motus_path, show_col_types = FALSE)
        rownames(df) <- df[[1]]
        df[[1]]      <- NULL

        # 5. optional alignment
        if (!raw) {
        aligned <- rename_and_align(
            counts_reprocessed = df,
            metadata          = metadata,
            scale             = scale,
            by_col            = "Sample_name",
            align             = align,
            study_name        = basename(local)
        )
        df <- aligned$reprocessed
        }

        # 6. numeric conversion + proportions
        df[]        <- lapply(df, as.numeric)
        proportions <- sweep(df, 1, rowSums(df), "/")

        # 7. simple taxonomy table from rownames
        tax_df <- tibble::tibble(taxa = rownames(df)) |>
        dplyr::mutate(taxa = stringr::str_trim(taxa)) |>
        tidyr::separate(
            taxa,
            into  = c("Kingdom","Phylum","Class","Order",
                    "Family","Genus","Species","Strain"),
            sep   = "\\s*;\\s*", extra = "drop", fill = "right"
        )
        rownames(tax_df) <- rownames(df)

        # 8. assign out
        mOTU3_counts       <- df
        mOTU3_proportions  <- proportions
        mOTU3_tax          <- tax_df
    }

    # 9. clean up only your private folder
    cleanup_tempfiles(temp_dir)
    }


    # ----- MetaPhlAn4 Reprocessed -----
    if (file.exists(metaphlan4_zip)) {
    # 1. private scratch folder
    temp_dir <- tempfile("mp4_")
    dir.create(temp_dir)

    # 2. locate the .tsv in the ZIP
    mp4_files    <- unzip(metaphlan4_zip, list = TRUE)
    mp4_filename <- mp4_files$Name[
                        grepl("\\.tsv$", mp4_files$Name, ignore.case = TRUE)
                    ][1]

    if (!is.na(mp4_filename)) {
        # 3. extract and capture full path
        unzipped  <- unzip(
                    metaphlan4_zip,
                    files     = mp4_filename,
                    exdir     = temp_dir,
                    overwrite = TRUE
                    )
        path      <- unzipped[1]

        # 4. read + set rownames
        df <- readr::read_tsv(path, show_col_types = FALSE)
        rownames(df) <- df[[1]]
        df[[1]]      <- NULL

        # 5. optional alignment
        if (!raw) {
        aligned <- rename_and_align(
            counts_reprocessed = df,
            metadata          = metadata,
            scale             = scale,
            by_col            = "Sample_name",
            align             = align,
            study_name        = basename(local)
        )
        df <- aligned$reprocessed
        }

        # 6. numeric + proportions
        df[]        <- lapply(df, as.numeric)
        proportions <- sweep(df, 1, rowSums(df), "/")

        # 7. taxonomy table
        tax_df <- tibble::tibble(taxa = rownames(df)) |>
        dplyr::mutate(taxa = stringr::str_trim(taxa)) |>
        tidyr::separate(
            taxa,
            into  = c("Kingdom","Phylum","Class","Order",
                    "Family","Genus","Species","Strain"),
            sep   = "\\s*;\\s*", extra = "drop", fill = "right"
        )
        rownames(tax_df) <- rownames(df)

        # 8. assign out
        MetaPhlAn4_counts      <- df
        MetaPhlAn4_proportions <- proportions
        MetaPhlAn4_tax         <- tax_df
    }

    # 9. tidy up
    cleanup_tempfiles(temp_dir)
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