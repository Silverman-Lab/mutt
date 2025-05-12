parse_2020_vieirasilva_nature_BMIS <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }
    if (!is.logical(raw) || length(raw) != 1) {
        stop("`raw` must be a single logical value (TRUE or FALSE)")
    }
    if (!is.logical(align) || length(align) != 1) {
        stop("`align` must be a single logical value (TRUE or FALSE)")
    }
    # Load needed libraries
    library(tidyverse)
    library(readxl)
    library(readr)

    # -------- local path -----------
    local                       <- file.path("2020_vieirasilva_nature_BMIS")

    # -------- file paths -----------
    original_proportions_zip    <- file.path(local, "qmp_genera_abundances_BMIS.tsv.zip")
    scale_zip                   <- file.path(local, "41586_2020_2269_MOESM3_ESM.xlsx.zip")
    metadata_zip                <- file.path(local, "SraRunTable (34).csv.zip")
    motus_zip                   <- file.path(local, "PRJEB37249_motus_merged.tsv.zip")
    metaphlan4_zip              <- file.path(local, "PRJEB37249_MetaPhlAn_merged.tsv.zip")

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

        # ----- Metadata -----
    metadata <- NA
    if (file.exists(metadata_zip)) {
        metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
        metadata_con <- unz(metadata_zip, metadata_csv)
        metadata <- read.csv(metadata_con, row.names = "Submitter_id") %>%
        as.data.frame() %>%
        rownames_to_column("Sample")
    }

    # ----- Scale ----- 
    scale <- NA
    if (file.exists(scale_zip)) {
        zinfo     <- unzip(scale_zip, list = TRUE)
        xlsx_name <- zinfo$Name[grepl("\\.xlsx$", zinfo$Name)][1]
        if (is.na(xlsx_name)) stop("No .xlsx file in ", scale_zip)
        temp_dir <- tempdir()
        unzip(scale_zip, files = xlsx_name, exdir = temp_dir)
        scale_path <- file.path(temp_dir, xlsx_name)
        COLUMN2 <- "Faecal microbial load (cells/g)"
        sheet2_full <- read_excel(scale_path, sheet = 3) %>% 
            as.data.frame()
        sheet2_full <- sheet2_full %>%
            rename(Sample = SampleID)
        sheet2_full <- sheet2_full %>%
            mutate(
            Sample = if_else(
                str_detect(Sample, "^x"),
                paste0("M0", Sample),
                Sample
            )
            )
        scale     <- sheet2_full %>% select(Sample, all_of(COLUMN2))
        metadata2 <- sheet2_full %>% select(-all_of(COLUMN2))
        cleanup_tempfiles(temp_dir)
    }

    scale = scale %>% 
        mutate(log2_FC_cells_g = ifelse(`Faecal microbial load (cells/g)` > 0, log2(`Faecal microbial load (cells/g)`), NA)) %>%
        mutate(log10_FC_cells_g = ifelse(`Faecal microbial load (cells/g)` > 0, log10(`Faecal microbial load (cells/g)`), NA))

    metadata <- metadata %>%
        full_join(metadata2, by = "Sample")

    # ----- original proportions --

    if (file.exists(original_proportions_zip)) {
        zinfo <- unzip(original_proportions_zip, list = TRUE)
        tsv_name <- zinfo$Name[grepl("\\.tsv$", zinfo$Name)][1]
        if (is.na(tsv_name)) stop("No .tsv file found inside ", original_proportions_zip)
        temp_dir <- tempdir()
        unzip(original_proportions_zip, files = tsv_name, exdir = temp_dir)
        tsv_path <- file.path(temp_dir, tsv_name)
        counts_original <- read.delim(tsv_path, row.names = 1, check.names = FALSE)
        counts_original <- as.data.frame(counts_original)
        counts_original <- counts_original %>%
            rename(Sample = SampleID)
        counts_original <- counts_original %>%
            mutate(
            Sample = if_else(
                str_detect(Sample, "^x"),
                paste0("M0", Sample),
                Sample
            )
        )
        if (!raw) {
            align = rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
            counts_original = align$counts_original
            original_names <- colnames(counts_original)
            counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
        }
        proportions_original <- sweep(counts_original, 1, rowSums(counts_original), FUN = "/")
        old_colnames <- colnames(proportions_original)
        tax_original <- data.frame(
        Taxa = old_colnames,
        stringsAsFactors = FALSE
        )
        cleanup_tempfiles(temp_dir)
    }

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
            by_col            = "SampleID",
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
            by_col            = "SampleID",
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