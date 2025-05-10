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
                aligned = rename_and_align(counts_reprocessed = df, metadata=metadata_meta_df, scale=scale_meta_df, by_col="Sample Name", align = align, study_name=basename(local))
                df = aligned$reprocessed
            }
            proportions <- sweep(df, 1, rowSums(df), FUN = "/")
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
                aligned = rename_and_align(counts_reprocessed = df, metadata=metadata_meta_df, scale=scale_meta_df, by_col="Sample Name", align = align, study_name=basename(local))
                df = aligned$reprocessed
            }
            proportions <- sweep(df, 1, rowSums(df), FUN = "/")
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

    # ----- Reprocessed counts from RDS ZIP -----
    if (file.exists(repro_counts_rds_zip)) {
    temp_rds <- tempfile(fileext = ".rds")
    unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

        rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
        if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(rds_files[1])) 

        # ----- Taxonomy reprocessed -----
        temp_tax <- tempfile(fileext = ".rds")
        unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

        tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
        if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa ----- 
        if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata_16s_df, scale=scale_16s_df, by_col="Sample Name", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
        }

        # taxa
        if (!raw) {
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
        }

        # proportions reprocessed
        proportions_reprocessed = counts_reprocessed
        proportions_reprocessed[-1] <- lapply(
            counts_reprocessed[-1],
            function(col) col / sum(col)
        )
    }

    # CAN DELETE LATER WHEN REPROCESSING FINISHES -- MY OLD REPROCESSING:
    counts_reprocessed    <- read_zipped_table(counts_16s_zip) 
    MetaPhlAn4_counts   <- read_zipped_table(counts_meta_zip) 
    if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata_16s_df, scale=scale_16s_df, by_col="Sample Name", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
        aligned = rename_and_align(counts_reprocessed = MetaPhlAn4_counts, metadata=metadata_meta_df, scale=scale_meta_df, by_col="Sample Name", align = align, study_name=basename(local))
        MetaPhlAn4_counts = aligned$reprocessed
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
    }

    return(list(
    counts = list(
                original = list(
                    amplicon = counts_original_16s,
                    shotgun = counts_original_meta
                ),
                reprocessed = list(
                    amplicon = counts_reprocessed,
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
                    amplicon = proportions_reprocessed,
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
                    amplicon = tax_reprocessed,
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