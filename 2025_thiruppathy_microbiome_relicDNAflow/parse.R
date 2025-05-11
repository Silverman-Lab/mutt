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
    metaphlan4_zip              <- file.path(local, "PRJNA1118035_MetaPhlAn_merged.tsv.zip")
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
        aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample", align = align)
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
                aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample", align = align)
                df = aligned$reprocessed
                original_names <- colnames(df)
                df <- as.data.frame(lapply(df, as.numeric), row.names = rownames(df), col.names = original_names, check.names = FALSE)
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
                aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample", align = align)
                df = aligned$reprocessed
                original_names <- colnames(df)
                df <- as.data.frame(lapply(df, as.numeric), row.names = rownames(df), col.names = original_names, check.names = FALSE)
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