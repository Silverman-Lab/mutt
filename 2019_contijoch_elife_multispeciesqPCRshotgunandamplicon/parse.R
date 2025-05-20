parse_2019_contijoch_elife_multispeciesqPCRshotgunandamplicon <- function(raw = FALSE, align = FALSE) {
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
    local <- file.path("2019_contijoch_elife_multispeciesqPCRshotgunandamplicon")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJNA413199_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA413199_dada2_taxa.rds.zip")
    motus_zip            <- file.path(local, "PRJNA413199_motus_merged.tsv.zip")
    metaphlan4_zip       <- file.path(local, "PRJNA413199_MetaPhlAn_merged.tsv.zip")
    scale_16s_zip        <- file.path(local, "Contijoch2019 16S_scale.csv.zip")
    scale_meta_zip       <- file.path(local, "Contijoch2019_scale.csv.zip")
    counts_16s_zip       <- file.path(local, "Contijoch_2019_16S.csv.zip")
    metadata_16s_zip     <- file.path(local, "Contijoch_2019_metadata.csv.zip")
    counts_meta_zip      <- file.path(local, "Contijoch_2019_shotgunmetagenomics.csv.zip")
    metadata_meta_zip    <- file.path(local, "Contijoch_metadata.csv.zip")

    # ----- Initialize everything as NA -----
    counts_original_16s <- NA
    proportions_original_16s <- NA
    tax_original_16s <- NA
    counts_original_meta <- NA
    proportions_original_meta <- NA
    tax_original_meta <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
    counts_reprocessed2 <- NA
    proportions_reprocessed2 <- NA
    tax_reprocessed2 <- NA
    mOTU3_counts <- NA
    mOTU3_proportions <- NA
    mOTU3_tax <- NA
    MetaPhlAn4_counts <- NA
    MetaPhlAn4_proportions <- NA
    MetaPhlAn4_tax <- NA

    # ------ original counts ------

    # ---- scale and metadata -----
    scale_16s_df     <- read_zipped_table(scale_16s_zip) %>% mutate(log2_Sample_microbial_density_ug_per_mg = ifelse(!is.na(Sample_microbial_density_ug_per_mg),
                                            Sample_microbial_density_ug_per_mg * log2(10),
                                            NA)) %>% rename(log10_Sample_microbial_density_ug_per_mg = Sample_microbial_density_ug_per_mg)
    scale_meta_df    <- read_zipped_table(scale_meta_zip) %>% mutate(log2_microbial_density = ifelse(!is.na(microbial_density),
                                            microbial_density * log2(10),
                                            NA)) %>% rename(log10_microbial_density = microbial_density)
    metadata_16s_df  <- read_zipped_table(metadata_16s_zip, row.names=NULL) %>% rename(Accession = Run)
    metadata_meta_df <- read_zipped_table(metadata_meta_zip, row.names=NULL) %>% rename(Accession = Run)

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
            scale             = scale_meta_df,
            by_col            = "Sample Name",
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
            scale             = scale_meta_df,
            by_col            = "Sample Name",
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
    if (file.exists(repro_counts_rds_zip)) {
        temp_rds <- tempfile("repro")
        dir.create(temp_rds)
        unzipped = unzip(repro_counts_rds_zip, exdir = temp_rds, overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(counts_file))
        cleanup_tempfiles(temp_rds)

        # ----- rdp16 -----
        if (!file.exists(file.path(local,"rdp16classified.csv.zip"))) {
        if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
            required_pkgs <- c("dada2", "Biostrings")
            missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
            if (length(missing_pkgs) > 0) {
                stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                    ". Please install them before running this function.")
            }
            seqs <- Biostrings::DNAStringSet(colnames(counts_reprocessed))
            rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
            tax_reprocessed2 = make_taxa_label(rdpclassified) 
            write.csv(tax_reprocessed2, file = file.path(local, "rdp16classified.csv"), row.names = TRUE)
            } else {
            stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
        }
        
        } else {
            tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
        }

        # ----- Taxonomy reprocessed -----
        unzipped = unzip(repro_tax_zip, exdir = temp_rds, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa ----- 
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata_16s_df, scale=scale_16s_df, by_col="Sample Name", align = align, study_name=basename(local))
            counts_reprocessed = aligned$reprocessed
            counts_reprocessed2 = aligned$reprocessed
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            matched_taxa2 <- tax_reprocessed2$Taxa[match(colnames(counts_reprocessed2), rownames(tax_reprocessed2))]
            colnames(counts_reprocessed) <- matched_taxa
            colnames(counts_reprocessed2) <- matched_taxa2
            counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
            counts_reprocessed2 <- collapse_duplicate_columns_exact(counts_reprocessed2)
            original_names <- colnames(counts_reprocessed)
            original_names2 <- colnames(counts_reprocessed2)
            counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
            counts_reprocessed2 <- as.data.frame(lapply(counts_reprocessed2, as.numeric), row.names = rownames(counts_reprocessed2), col.names = original_names2, check.names = FALSE)
            proportions_reprocessed2 <- sweep(counts_reprocessed2, 1, rowSums(counts_reprocessed2), '/')
        }
        # proportions reprocessed
        proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
        cleanup_tempfiles(temp_rds)
    }

    # CAN DELETE LATER WHEN REPROCESSING FINISHES -- MY OLD REPROCESSING:
    counts_reprocessed    <- read_zipped_table(counts_16s_zip) 
    MetaPhlAn4_counts     <- read_zipped_table(counts_meta_zip) 
    MetaPhlAn4_tax        <- data.frame(Taxa = colnames(MetaPhlAn4_counts))
    tax_reprocessed       <- data.frame(Taxa = colnames(counts_reprocessed))
    if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata_16s_df, scale=scale_16s_df, by_col="Sample Name", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
        original_names <- colnames(counts_reprocessed)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
        aligned = rename_and_align(counts_reprocessed = MetaPhlAn4_counts, metadata=metadata_meta_df, scale=scale_meta_df, by_col="Sample Name", align = align, study_name=basename(local))
        MetaPhlAn4_counts = aligned$reprocessed
        original_names <- colnames(MetaPhlAn4_counts)
        MetaPhlAn4_counts <- as.data.frame(lapply(MetaPhlAn4_counts, as.numeric), row.names = rownames(MetaPhlAn4_counts), col.names = original_names, check.names = FALSE)
    }
    proportions_reprocessed <- sweep(counts_reprocessed, MARGIN = 1,STATS  = rowSums(counts_reprocessed), FUN = "/")
    MetaPhlAn4_proportions <- sweep(MetaPhlAn4_counts, MARGIN = 1,STATS  = rowSums(MetaPhlAn4_counts), FUN = "/")


    if (!raw) {
        counts_original_16s = fill_na_zero_numeric(counts_original_16s)
        counts_original_meta = fill_na_zero_numeric(counts_original_meta)
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_original_16s = fill_na_zero_numeric(proportions_original_16s)
        proportions_original_meta = fill_na_zero_numeric(proportions_original_meta)
        proportions = fill_na_zero_numeric(proportions)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
        MetaPhlAn4_counts = fill_na_zero_numeric(MetaPhlAn4_counts)
        mOTU3_proportions = fill_na_zero_numeric(mOTU3_proportions)
        MetaPhlAn4_proportions = fill_na_zero_numeric(MetaPhlAn4_proportions)
        counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
        proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    }

    return(list(
    counts = list(
                original = list(
                    amplicon = counts_original_16s,
                    shotgun = counts_original_meta
                ),
                reprocessed = list(
                    amplicon = list(
                        rdp19 = counts_reprocessed,
                        rdp16 = counts_reprocessed2
                    ),
                    shotgun = list(
                        mOTU3 = mOTU3_counts,
                        MetaPhlan4 = MetaPhlAn4_counts
                    )
                )

    ),
    proportions = list(
                original = list(
                    amplicon = proportions_original_16s,
                    shotgun = proportions_original_meta
                ),
                reprocessed = list(
                    amplicon = list(
                        rdp19 = proportions_reprocessed,
                        rdp16 = proportions_reprocessed2
                    ),
                    shotgun = list(
                        mOTU3 = mOTU3_proportions,
                        MetaPhlan4 = MetaPhlAn4_proportions
                    )
                )
    ),
    tax = list(
                original = list(
                    amplicon = tax_original_16s,
                    shotgun = tax_original_meta
                ),
                reprocessed = list(
                    amplicon = list(
                        rdp19 = tax_reprocessed,
                        rdp16 = tax_reprocessed2
                    ),
                    shotgun = list(
                        mOTU3 = mOTU3_tax,
                        MetaPhlan4 = MetaPhlAn4_tax
                    )
                )
    ),
    scale = list(
                amplicon = scale_16s_df,
                shotgun = scale_meta_df
    ),
    metadata = list(
                amplicon = metadata_16s_df,
                shotgun = metadata_meta_df
    )
    ))
}