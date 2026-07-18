parse_2022_zaramela_msystems_synDNA <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }
    if (!is.logical(raw)) {
        stop("raw must be a boolean")
    }
    if (!is.logical(align)) {
        stop("align must be a boolean")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2022_zaramela_msystems_synDNA")

    # ---- initialize ----
    mOTU3_counts = NA
    mOTU3_proportions = NA
    mOTU3_tax = NA
    MetaPhlAn4_counts = NA
    MetaPhlAn4_proportions = NA
    MetaPhlAn4_tax = NA
    counts = NA
    proportions = NA
    tax = NA

    # ----- File paths -----
    motus_zip            <- file.path(local, "PRJNA870099_motus_merged.tsv.zip")
    metaphlan4_zip       <- file.path(local, "PRJNA870099_MetaPhlAn_merged_counts.tsv.zip")
    metadata_SRA_zip     <- file.path(local, "SraRunTable (34).csv.zip")
    metadatareadin       <- file.path(local, "synDNA_metadata_updated.tsv.zip")
    countsreadin         <- file.path(local, "HMP_frequency_table.tsv.zip")
    scalereadin          <- file.path(local, "Saliva_Flow_data.tsv.zip")
    taxreadin            <- file.path(local, "taxaID_length_fulllineage.tsv.zip")

    # ----- scale and metadata -----
    metadata1 <- read_zipped_table(metadata_SRA_zip, row.names=NULL) %>%
                    rename(Accession = Run, Sample = `Sample Name`)
    metadata2 <- read_zipped_table(metadatareadin, sep="\t")
    metadata2 <- metadata2 %>% rename(SampleID = Sample) %>% mutate(Sample = paste(SampleID, Pool, sep = "_"))
    metadata = full_join(metadata1, metadata2, by="Sample")
    scale <- read_zipped_table(scalereadin, row.names=NULL, sep="\t") %>%
             rename(Sample = SampleID) 

    metadata = left_join(metadata, scale %>% select(c("Sample", "gender", "saliva_weight_g", 
                "saliva_volume_mL_collected_in_5_min", "saliva_flow_rate_mL_per_min", "FC_avg_cells_per_ul", "FC_avg_cells_5_min")), by=c("SampleID" = "Sample"))
    scale = metadata %>% select(ID, Sample, FC_avg_cells_per_ul, FC_avg_cells_5_min)  %>% 
        mutate(log2_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log2(FC_avg_cells_per_ul), NA)) %>%
        mutate(log2_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log2(FC_avg_cells_5_min), NA)) %>%
        mutate(log10_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log10(FC_avg_cells_per_ul), NA)) %>%
        mutate(log10_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log10(FC_avg_cells_5_min), NA))

    scale = scale %>% filter(!is.na(FC_avg_cells_per_ul)) %>% select(-Sample)
    metadata = metadata %>% select(-c("FC_avg_cells_per_ul", "FC_avg_cells_5_min"))

    metadata <- metadata %>%
        mutate(across(c(gender, Pool, Strand), factor))

    # ----- counts and proportions and taxa -----
    counts <- read_zipped_table(countsreadin, row.names=NULL, sep="\t") %>%
                group_by(OTUID) %>%  
                summarise(across(where(is.numeric), sum, na.rm = TRUE))
    counts <- as.data.frame(counts)
    rownames(counts) <- counts$OTUID
    counts$OTUID <- NULL
    counts = as.data.frame(t(counts))

    tmp_dir <- tempfile("tax")
    dir.create(tmp_dir)
    tmp_file <- unzip(taxreadin, exdir = tmp_dir)[1]
    tax <- read_tsv(tmp_file, quote = "", col_types = cols(.default = "c")) %>%
    rename(Kingdom = superkingdom, Phylum = phylum, Class = class, 
            Order = order, Family = family, Genus = genus, Species = species, 
            OTUID = GenomeID) %>%
    mutate(
        MetaPhlAn4_lineage = paste(
            paste0("k__", Kingdom),
            paste0("p__", Phylum),
            paste0("c__", Class),
            paste0("o__", Order),
            paste0("f__", Family),
            paste0("g__", Genus),
            paste0("s__", Species),
            sep = "|"
        )
    )

    tax = as.data.frame(make_taxa_label(tax))
    rownames(tax) = tax$OTUID
    cleanup_tempfiles(tmp_dir)

    if (!raw) {
        aligned = rename_and_align(counts_original = counts, metadata = metadata, scale = scale, by_col = "ID", align = align, study_name = basename(local))
        counts = aligned$counts_original
        matched_taxa <- tax$MetaPhlAn4_lineage[match(colnames(counts), rownames(tax))]
        colnames(counts) <- matched_taxa
        counts <- collapse_duplicate_columns_exact(counts)
        original_names <- colnames(counts)
        counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
    }

    # --- Compute proportions from counts ---
    proportions <- sweep(counts, MARGIN = 1,STATS  = rowSums(counts), FUN = "/")

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
            by_col            = "ID",
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
            unzipped <- unzip(
                metaphlan4_zip,
                files     = mp4_filename,
                exdir     = temp_dir,
                overwrite = TRUE
            )
            path <- unzipped[1]

            # 4. read + set rownames
            df <- readr::read_tsv(path, show_col_types = FALSE) %>%
                as.data.frame() %>%
                column_to_rownames("clade") %>%
                t() %>%
                as.data.frame()

            # 5. optional alignment
            if (!raw) {
                aligned <- rename_and_align(
                    counts_reprocessed = df,
                    metadata           = metadata,
                    scale              = scale,
                    by_col             = "ID",
                    align              = align,
                    study_name         = basename(local)
                )
                df <- aligned$reprocessed
            }

            # 6. numeric conversion
            df[] <- lapply(df, as.numeric)


            tax_df <- data.frame(taxa = colnames(df)) %>%
                mutate(taxa = str_trim(taxa)) %>%
                separate(
                    taxa,
                    into  = c("Kingdom", "Phylum", "Class", "Order",
                            "Family", "Genus", "Species", "Strain"),
                    sep   = "\\|",
                    extra = "drop",
                    fill  = "right"
                ) %>%
                # remove the leading “letter__” prefix (e.g. “k__”, “p__”…)
                mutate(across(Kingdom:Strain, ~ str_remove(.x, "^[a-z]__")))

            rownames(tax_df) <- colnames(df)
            tax_df <- make_taxa_label(tax_df)

            # 7. keep only columns that contain "s__" but do NOT contain "t__"
            df <- df[, grepl("s__", colnames(df)) & !grepl("t__", colnames(df)), drop = FALSE]


            # 8. compute proportions
            proportions <- sweep(df, 1, rowSums(df), "/")

            # 10. assign outputs
            MetaPhlAn4_counts      <- df
            MetaPhlAn4_proportions <- proportions
            MetaPhlAn4_tax         <- tax_df
        }

        # 11. tidy up
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

    # ----- Return -----
    return(list(
        counts = list(
        original = counts,
        reprocessed = list(
            mOTU3 = mOTU3_counts,
            MetaPhlAn4 = MetaPhlAn4_counts
        )
        ),
        proportions = list(
        original = proportions,
        reprocessed = list(
            mOTU3 = mOTU3_proportions,
            MetaPhlAn4 = MetaPhlAn4_proportions
        )
        ),
        tax = list(
        original = tax,
        reprocessed = list(
            mOTU3 = mOTU3_tax,
            MetaPhlAn4 = MetaPhlAn4_tax
        )
        ),
        scale = scale,
        metadata = metadata
    ))
}
