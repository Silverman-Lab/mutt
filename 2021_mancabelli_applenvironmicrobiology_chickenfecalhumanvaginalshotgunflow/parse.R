parse_2021_mancabelli_applenvironmicrobiology_chickenfecalhumanvaginalshotgunflow <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }
    if (!is.logical(raw)) {
        stop("raw must be a logical value")
    }
    if (!is.logical(align)) {
        stop("align must be a logical value")
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
    metadata_df  <- read_zipped_table(metadata_zip) %>% rename(Accession = Run)

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
            metadata          = metadata_df,
            scale             = scale_df,
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
            metadata          = metadata_df,
            scale             = scale_df,
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
        counts_original = fill_na_zero_numeric(counts_original)
        mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
        proportions_original = fill_na_zero_numeric(proportions_original)
        MetaPhlAn4_counts = fill_na_zero_numeric(MetaPhlAn4_counts)
        mOTU3_proportions = fill_na_zero_numeric(mOTU3_proportions)
        MetaPhlAn4_proportions = fill_na_zero_numeric(MetaPhlAn4_proportions)
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