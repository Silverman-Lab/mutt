parse_2025_thiruppathy_microbiome_relicDNAflow <- function(raw = FALSE, align = FALSE) {
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
    local                       <- file.path("2025_thiruppathy_microbiome_relicDNAflow")

    # -------- file paths -----------
    original_counts_zip         <- file.path(local, "Final_rpca_all.tsv.zip")
    scale_zip                   <- file.path(local, "scale.csv.zip")
    metadata_three_zip          <- file.path(local, "metadata.csv.zip")
    metadata_two_zip            <- file.path(local, "40168_2025_2063_MOESM3_ESM.csv.zip")
    metadata_zip                <- file.path(local, "SraRunTable (38).csv.zip")
    motus_zip                   <- file.path(local, "PRJNA1118035_motus_merged.tsv.zip")
    metaphlan4_zip              <- file.path(local, "PRJNA1118035_MetaPhlAn_merged_counts.tsv.zip")
    tax_zip                     <- file.path(local, "SMGC_60perc_covered_genomes.txt.zip")

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
    sra         = read_zipped_table(metadata_zip, row.names=NULL) %>% as.data.frame() %>% rename(Accession = Run, Sample = `Sample Name`)
    metadata1   = read_zipped_table(metadata_three_zip, row.names=NULL) %>% as.data.frame()
    metadata2   = read_zipped_table(metadata_two_zip, row.names=NULL) %>% as.data.frame()
    metadata    = sra %>% full_join(metadata1, by = "Sample") %>% full_join(metadata2, by = "Sample")

    # ----- Scale -----
    scale       = read_zipped_table(scale_zip, row.names=NULL) %>% as.data.frame() %>% mutate(log2_Total_Cells_Per_sq_cm = ifelse(Total_Cells_Per_sq_cm > 0,
                                    log2(Total_Cells_Per_sq_cm),NA)) %>% mutate(log10_Total_Cells_Per_sq_cm = ifelse(Total_Cells_Per_sq_cm > 0, 
                                    log10(Total_Cells_Per_sq_cm),NA))
    colnames(scale) <- gsub("\\n", "", colnames(scale))
    mergedmeta  = scale %>% select(c("Sample", "Dilution", "Bead_Number", "Bacteria_Single_Cell_Number", "AccuCount_Particle_Number_30uL", "Volume_Run_In_Machine_uL", "Swab_Area_sq_cm", "PBS_Volume_For_Swab_uL"))                                          
    scale       = scale %>% select(-c( "Dilution", "Bead_Number", "Bacteria_Single_Cell_Number","AccuCount_Particle_Number_30uL", "Volume_Run_In_Machine_uL", "Swab_Area_sq_cm", "PBS_Volume_For_Swab_uL"))  
    metadata    = metadata %>% full_join(mergedmeta, by = "Sample")

    # ----- original counts, tax, proportions -----
    counts_original = read_zipped_table(original_counts_zip, sep = "\t", row.names=NULL) %>% 
                        rename(OTU_ID = `OTU ID`) %>% t() %>% as.data.frame() 
    colnames(counts_original) <- as.character(counts_original[1, ]) 
    counts_original <- counts_original[-1, ]  

    # --- taxa ---
    tax_original = read_zipped_table(tax_zip, sep = "\t", row.names=NULL) %>% as.data.frame() %>% 
                    rename(OTU_ID = gotu, taxonomy = strain) %>% select(-c("covered_length", "total_length"))
    tax_original$ogtaxonomy = tax_original$taxonomy                
    tax_original <- tax_original %>% separate(taxonomy,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
           sep = ";", fill = "right") %>% mutate(across(Kingdom:Strain, ~ gsub("^[a-z]__+", "", .)))                
    tax_original = make_taxa_label(tax_original)
    rownames(tax_original) <- tax_original$OTU_ID
    
    if (!raw) {
        aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
        counts_original = aligned$counts_original
        matched_taxa <- tax_original$ogtaxonomy[match(colnames(counts_original), rownames(tax_original))]
        colnames(counts_original) <- matched_taxa
        counts_original <- collapse_duplicate_columns_exact(counts_original)
        original_names <- colnames(counts_original)
        counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
       
    }

    # ------ proportions from counts ------
    proportions_original <- sweep(counts_original, MARGIN = 1,STATS  = rowSums(counts_original), FUN = "/")

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
            by_col            = "Sample",
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
                    by_col             = "Sample",
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