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
                aligned = rename_and_align(counts_original = df, metadata = metadata, scale = scale, by_col = "Sample_ID", align = align, study_name = basename(local))
                df = aligned$counts_original
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
                aligned = rename_and_align(counts_original = df, metadata = metadata, scale = scale, by_col = "Sample_ID", align = align, study_name = basename(local))
                df = aligned$counts_original
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
    if (all(file.exists(repro_counts_zips))) {
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
            aligned <- rename_and_align(counts_reprocessed=counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample_ID", align = align, study_name = basename(local))
            counts_reprocessed <- aligned$reprocessed
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