parse_2020_barlow_naturecommunications_miceGI <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
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

    # Load libraries
    library(tidyverse)
    library(readxl)
    library(stringr)
    library(readr)

    # ----- Local base directory -----
    local <- file.path("2020_barlow_naturecommunications_miceGI")

    # ----- File paths -----
    sra_metadata_zip     <- file.path(local, "SraRunTable (32).csv.zip")
    counts               <- NA  # No original counts provided
    sraplusload_zip      <- file.path(local, "SraRunTable_total_loads.csv.zip")
    proportions_zip      <- file.path(local, "Relative_Abundance_Table.csv.zip")
    scale_zip            <- file.path(local, "absolute_abundance_calculated.csv.zip")
    repro_counts_rds_zip <- file.path(local, "PRJNA575097_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA575097_dada2_taxa.rds.zip")

    counts_reprocessed2 <- NA
    tax_reprocessed2 <- NA
    proportions_reprocessed2 <- NA
    counts_reprocessed <- NA
    tax_reprocessed <- NA
    proportions_reprocessed <- NA
    counts = NA
    proportions = NA    
    tax_original = NA

    # ----- metadata ------
    sra_metadata = read_zipped_table(sra_metadata_zip, row.names = NULL) %>% 
                    mutate(Sample = str_replace_all(`Sample Name`, "_", " ")) %>% 
                    rename(Accession = Run)
    
    # --- Scale (ddPCR) ---
    scale_csv <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_csv, exdir = tempdir(), overwrite = TRUE)
    scale <- read.csv(scale_path, row.names = NULL, stringsAsFactors = FALSE)
    metadata <- full_join(sra_metadata, scale, by = "Sample") 
    metadata = read_zipped_table(sraplusload_zip, row.names = NULL) %>% rename(Accession = Run) %>%
                    full_join(metadata, by = "Accession")
    dups <- duplicated(as.list(metadata))
    metadata <- metadata[, !dups] %>% 
        as.data.frame() 
    scale = metadata %>% select(c("Sample" ,"Accession", "Total Load (16S Copies/g)")) %>%
            mutate(log2_total_16S_copies_g = ifelse(`Total Load (16S Copies/g)` > 0, log2(`Total Load (16S Copies/g)`), NA)) %>% 
            mutate(log10_total_16S_copies_g = ifelse(`Total Load (16S Copies/g)` > 0, log10(`Total Load (16S Copies/g)`), NA))
    metadata = metadata %>% select(-c("Total Load (16S Copies/g)"))

        # --- proportions ---
    prop_csv <- unzip(proportions_zip, list = TRUE)$Name[1]
    prop_path <- unzip(proportions_zip, files = prop_csv, exdir = tempdir(), overwrite = TRUE)
    proportions <- read.csv(prop_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

    # metadatabackup <- proportions %>%
    #     rownames_to_column("Sample") %>%
    #     select(c(Sample, Diet, Site, Day, mouse, Cage)) %>% 
    #     mutate(Sample = str_replace_all(Sample, "_", " "))

    proportions <- proportions %>% select(-c(Diet, Site, Day, mouse, Cage))

    if (!raw) {
        aligned = rename_and_align(proportions_original = proportions, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
        proportions = aligned$proportions_original
        original_names <- colnames(proportions)
        proportions <- as.data.frame(lapply(proportions, as.numeric), row.names = rownames(proportions), col.names = original_names, check.names = FALSE)
    }
    original_taxa <- colnames(proportions)

    # --- Create taxa dataframe (lookup table) ---
    tax <- data.frame(taxonomy = original_taxa, stringsAsFactors = FALSE)
    rownames(tax) <- tax$taxonomy
    tax$ogtaxonomy <- tax$taxonomy
    tax_original <- tax %>%
        separate(taxonomy,
                    into = paste0("D", 0:6),
                    sep = ";", fill = "right", remove = FALSE) %>%
        mutate(across(starts_with("D"), ~ gsub("^D_\\d+__", "", .))) %>%
        mutate(across(starts_with("D"), ~ ifelse(. == "" | . == "__" | is.na(.), "unclassified", .))) %>%
        transmute(
            Kingdom = D0,
            Phylum  = D1,
            Class   = D2,
            Order   = D3,
            Family  = D4,
            Genus   = D5,
            Species = D6
        )

    tax_original = make_taxa_label(tax_original) %>% as.data.frame(stringsAsFactors = FALSE)
    proportions <- as.data.frame(proportions, stringsAsFactors = FALSE)
    if (!raw) {
        matched_taxa <- tax_original$Taxa[match(colnames(proportions), rownames(tax_original))]
        colnames(proportions) <- matched_taxa
        proportions <- collapse_duplicate_columns_exact(proportions)
        original_names <- colnames(proportions)
        proportions <- as.data.frame(lapply(proportions, as.numeric), row.names = rownames(proportions), col.names = original_names, check.names = FALSE)
    }
    
    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
        temp_dir <- file.path(tempdir(), "repro")
        dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
        unzipped = unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(counts_file))

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
        unzipped = unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))

        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
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
        cleanup_tempfiles(temp_dir)
    }

    if (!raw) {
      counts = fill_na_zero_numeric(counts)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions = fill_na_zero_numeric(proportions)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
      proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    }

    # ---- return ----
    return(list(
        counts = list(
            original = counts,
            reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
        ),
        tax = list(
            original = tax_original,
            reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
        ),
        proportions = list(
            original = proportions,
            reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
        ),
        metadata = metadata,
        scale = scale
    ))
}
