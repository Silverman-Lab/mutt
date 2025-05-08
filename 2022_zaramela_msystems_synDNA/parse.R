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

    # ----- File paths -----
    motus_zip            <- file.path(local, "PRJNA940499_motus_merged.tsv.zip")
    metaphlan4_zip       <- file.path(local, "PRJNA940499_MetaPhlAn_merged.tsv.zip")
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
                "saliva_volume_mL_collected_in_5_min", "saliva_flow_rate_mL_per_min")), by="Sample")
    scale = scale %>% select(-c("gender"))  %>% 
        mutate(log2_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log2(FC_avg_cells_per_ul), NA)) %>%
        mutate(log2_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log2(FC_avg_cells_5_min), NA)) %>%
        mutate(log10_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log10(FC_avg_cells_per_ul), NA)) %>%
        mutate(log10_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log10(FC_avg_cells_5_min), NA))


    # ----- counts and proportions and taxa -----

    counts <- read_zipped_table(countsreadin, row.names=NULL, sep="\t") %>%
                group_by(OTUID) %>%  
                summarise(across(where(is.numeric), sum, na.rm = TRUE))
    counts <- as.data.frame(counts)
    rownames(counts) <- counts$OTUID
    counts$OTUID <- NULL
    counts = as.data.frame(t(counts))

    tmp_file <- unzip(taxreadin, exdir = tempdir())[1]
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

    tax = make_taxa_label(tax)
    rownames(tax) = tax$OTUID

    if (!raw) {
        matched_taxa <- tax$MetaPhlAn4_lineage[match(colnames(counts), rownames(tax))]
        colnames(counts) <- matched_taxa
        counts <- as.data.frame(t(rowsum(t(counts), group = colnames(counts))))
    }

    # --- Compute proportions from counts ---
    proportions <- sweep(counts, MARGIN = 1,STATS  = rowSums(counts), FUN = "/")

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
                aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
                df = aligned$reprocessed
            }
            proportions <- apply(df, 2, function(col) col / sum(col))
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
                aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
                df = aligned$reprocessed
            }
            proportions <- apply(df, 2, function(col) col / sum(col))
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
