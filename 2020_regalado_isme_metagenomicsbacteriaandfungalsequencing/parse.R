parse_2020_regalado_isme_metagenomicsbacteriaandfungalsequencing <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }
    if (!is.logical(raw) || length(raw) != 1) {
        stop("`raw` must be a single logical value (TRUE or FALSE)")
    }
    if (!is.logical(align) || length(align) != 1) {
        stop("`align` must be a single logical value (TRUE or FALSE)")
    }

    library(stringr)
    library(tidyverse)
    # ----- Local base directory -----
    local <- file.path("2020_regalado_isme_metagenomicsbacteriaandfungalsequencing")

    # ----- File paths -----
    motus_zip                   <- file.path(local, "PRJNA31530_motus_merged.tsv.zip")
    metaphlan4_zip              <- file.path(local, "PRJNA31530_MetaPhlAn_merged.tsv.zip")
    metadata_zip                <- file.path(local, "SraRunTable (38).csv.zip")
    ITS_metadata_zip            <- file.path(local, "ITS_blackblue_metadata_v2.txt.zip")
    ITS_og_otu_zip              <- file.path(local, "AgITS1blackblue_READ1_271_all_Zotutab.txt.zip")
    ITS_og_tax_zip              <- file.path(local, "AgITS1blackblue_READ1_271_all_Zotus.tax.zip")
    16S_og_otu_zip              <- file.path(local, "515_blackblue_all_Zotutab_20181206.txt.zip")
    16S_og_tax_zip              <- file.path(local, "515_blackblue_all_zotus.tax.zip")
    16S_og_metadata_zip         <- file.path(local, "515_blackblue_metadata.txt.zip")
    metagenomic_og_bacteria_zip <- file.path(local, "BacteriaGenusRaw.txt.zip")
    metagenomic_og_fungi_zip    <- file.path(local, "FungiGenusRaw.txt.zip")
    scale_16s_zip               <- file.path(local, "total_seq.txt.zip")
    metadataconnection          <- file.path(local, "Metadata_PRJEB31530.csv.zip")

    # ----- Initialize everything as NA -----
    counts_original_16s <- NA
    proportions_original_16s <- NA
    tax_original_16s <- NA
    counts_original_ITS <- NA
    proportions_original_ITS <- NA
    tax_original_ITS <- NA
    counts_original_meta <- NA
    proportions_original_meta <- NA
    tax_original_meta <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
    mOTU3_counts <- NA
    mOTU3_proportions <- NA
    mOTU3_tax <- NA
    MetaPhlAn4_counts <- NA
    MetaPhlAn4_proportions <- NA
    MetaPhlAn4_tax <- NA
    counts_ITS_reprocessed <- NA
    proportions_ITS_reprocessed <- NA
    tax_ITS_reprocessed <- NA


    repro_counts_zips <- c(
    file.path(local, "PRJNA31530_dada2_counts.rds.zip"),
    file.path(local, "PRJNA31530_SILVA_counts.rds.zip")
    )

    repro_tax_zips <- c(
    file.path(local, "PRJNA31530_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA31530_SILVA_taxa.rds.zip")
    )

    # ---- scale and metadata -----
    scale_16s_df        <- read_zipped_table(scale_16s_zip, sep = "\t", row.names = NULL) %>% mutate(log2_qPCR = ifelse(TotalSeq > 0, log2(TotalSeq), NA)) %>% 
                                                             mutate(log10_qPCR = ifelse(TotalSeq > 0, log10(TotalSeq), NA)) 
    metadata_16s_df     <- read_zipped_table(metadata_16s_zip)
    metadata_meta_df    <- read_zipped_table(metadata_meta_zip)
    metadata_ITS_df     <- read_zipped_table(metadata_ITS_zip)
    metadata_connection <- read_zipped_table(metadataconnection)

    # ------ original counts ------
    counts_original_16s <- read_zipped_table(16S_og_otu_zip)
    counts_original_ITS <- read_zipped_table(ITS_og_otu_zip)
    counts_original_meta <- read_zipped_table(metagenomic_og_bacteria_zip)
    counts_original_meta <- read_zipped_table(metagenomic_og_fungi_zip)

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
            metadata          = metadata_meta_df,
            scale             = scale_16s_df,
            by_col            = "Sample_ID",
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
            metadata          = metadata_meta_df,
            scale             = scale_16s_df,
            by_col            = "Sample_ID",
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


    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(repro_counts_zips))) {
        temp_dir <- tempdir("repro")
        dir.create(temp_dir)
        unzipped = unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(counts_file))

        # ----- Taxonomy reprocessed -----
        unzipped = unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))

        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            aligned <- rename_and_align(counts_reprocessed=counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample_ID", align = align, study_name = basename(local))
            counts_reprocessed <- aligned$reprocessed
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed = collapse_duplicate_columns_exact(counts_reprocessed)
            original_names <- colnames(counts_reprocessed)
            counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
        }

        # proportions reprocessed
        proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
        cleanup_tempfiles(temp_dir)
    }


    if (!raw) {
        counts = fill_na_zero_numeric(counts)
        mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
        proportions = fill_na_zero_numeric(proportions)
        MetaPhlAn4_counts = fill_na_zero_numeric(MetaPhlAn4_counts)
        mOTU3_proportions = fill_na_zero_numeric(mOTU3_proportions)
        MetaPhlAn4_proportions = fill_na_zero_numeric(MetaPhlAn4_proportions)
        counts_ITS_reprocessed = fill_na_zero_numeric(counts_ITS_reprocessed)
        proportions_ITS_reprocessed = fill_na_zero_numeric(proportions_ITS_reprocessed)
    }

    return(list(
        counts = list(
                    original = counts,
                    reprocessed = list(
                        amplicon = counts_reprocessed,
                        shotgun = list(
                            mOTU3 = mOTU3_counts,
                            MetaPhlan4 = MetaPhlAn4_counts
                        ),
                        ITS = counts_ITS_reprocessed
                    )

        ),
        proportions = list(
                    original = proportions,
                    reprocessed = list(
                        amplicon = proportions_reprocessed,
                        shotgun = list(
                            mOTU3 = mOTU3_proportions,
                            MetaPhlan4 = MetaPhlAn4_proportions
                        ),
                        ITS = proportions_ITS_reprocessed
                    )
        ),
        tax = list(
                    original = tax,
                    reprocessed = list(
                        amplicon = tax_reprocessed,
                        shotgun = list(
                            mOTU3 = mOTU3_tax,
                            MetaPhlan4 = MetaPhlAn4_tax
                        ),
                        ITS = tax_ITS_reprocessed
                    )

        ),
        scale = list(
                    amplicon = scale_16s_df,
                    shotgun = scale_meta_df,
                    ITS = scale_ITS_df
        ),
        metadata = list(
                    amplicon = metadata_16s_df,
                    shotgun = metadata_meta_df,
                    ITS = metadata_ITS_df
        )
    ))
}