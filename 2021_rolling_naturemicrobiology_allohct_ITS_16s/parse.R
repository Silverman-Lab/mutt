parse_2021_rolling_naturemicrobiology_allohct_ITS_16s <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }
    if (!is.logical(raw)) {
        stop("raw must be a logical value")
    }
    if (!is.logical(align)) {
        stop("align must be a logical value")
    }
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2021_rolling_naturemicrobiology_allohct_ITS_16s")

    # ----- File paths -----
    scale_zip     <- file.path(local, "scale.csv.zip")
    metadata_zip  <- file.path(local, "metadata.csv.zip")
    sra_zips       <- c(file.path(local, "SraRunTable (43).csv.zip"),
                       file.path(local, "SraRunTable (44).csv.zip"),
                       file.path(local, "SraRunTable (45).csv.zip"),
                       file.path(local, "SraRunTable (46).csv.zip"),
                       file.path(local, "SraRunTable (47).csv.zip"),
                       file.path(local, "SraRunTable (48).csv.zip"),
                       file.path(local, "SraRunTable (49).csv.zip"),
                       file.path(local, "SraRunTable (50).csv.zip"),
                       file.path(local, "SraRunTable (51).csv.zip"),
                       file.path(local, "SraRunTable (52).csv.zip"),
                       file.path(local, "SraRunTable (53).csv.zip"),
                       file.path(local, "SraRunTable (54).csv.zip"),
                       file.path(local, "SraRunTable (55).csv.zip")
                      )
    culture_cfu_zip <- file.path(local, "culture_cfu.csv.zip")
    repro_counts_zips <- c(file.path(local, "PRJNA746305_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545312_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545313_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545314_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545315_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545316_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545317_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545318_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545319_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA545320_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA579121_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA548153_dada2_counts.rds.zip"),
                           file.path(local, "PRJNA607574_dada2_counts.rds.zip")
    )

    repro_tax_zips <- c(file.path(local, "PRJNA746305_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545312_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545313_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545314_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545315_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545316_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545317_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545318_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545319_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA545320_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA579121_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA548153_dada2_taxa.rds.zip"),
                        file.path(local, "PRJNA607574_dada2_taxa.rds.zip")
    )


    # ---- initialize ----
    scale <- NULL
    metadata <- NULL
    counts <- NULL
    tax <- NULL
    proportions <- NULL
    counts_reprocessed <- NULL
    tax_reprocessed <- NULL
    proportions_reprocessed <- NULL
    culture_cfu <- NULL
    sra <- NULL
    counts_reprocessed2 <- NULL
    tax_reprocessed2 <- NULL
    proportions_reprocessed2 <- NULL
    combined_counts <- NULL
    combined_props <- NULL
    combined_taxa <- NULL
    combined_counts2 <- NULL
    combined_props2 <- NULL
    combined_taxa2 <- NULL

    # ---- scale ----
    scale <- read_zipped_table(scale_zip, row.names=NULL)

    # ---- culture cfu ----
    culture_cfu <- read_zipped_table(culture_cfu_zip, row.names=NULL)

    # ---- metadata ----
    metadata <- read_zipped_table(metadata_zip, row.names=NULL)

    # ---- sra ----
    sra_tables <- lapply(sra_zips, function(zip) {
        if (file.exists(zip)) {
        read_zipped_table(zip, row.names = NULL) %>% rename(Accession = Run)
        } else {
        NULL
        }
    })
    
    # Combine all tables keeping all columns
    sra <- bind_rows(sra_tables[!sapply(sra_tables, is.null)])


    # ---- counts, proportions, tax ----   
    
    # ---- reprocessed counts, proportions, tax ----
    all_counts = NA
    all_props = NA
    all_taxa = NA
    all_counts2 = NA
    all_props2 = NA
    all_taxa2 = NA
    all_counts <- list()
    all_props  <- list()
    all_taxa   <- list()
    all_counts2 <- list()
    all_props2 <- list()
    all_taxa2 <- list()

    if (all(file.exists(repro_counts_zips)) && all(file.exists(repro_tax_zips))) {
        for (i in seq_along(repro_counts_zips)) {
        counts_zip <- repro_counts_zips[i]
        tax_zip    <- repro_tax_zips[i]

        # ----- Study prefix -----
        study_prefix <- gsub("_dada2_counts\\.rds\\.zip$", "", basename(counts_zip))

        # ----- Unzip and read counts -----
        temp_dir <- tempfile("repro")
        dir.create(temp_dir)
        unzipped = unzip(counts_zip, exdir = temp_dir, overwrite = TRUE)
        counts_files = unzipped[grepl("_counts\\.rds$", unzipped)][1]
        if (length(counts_files) == 0) stop("No *_counts.rds file found in: ", counts_zip)
        counts_reprocessed <- as.data.frame(readRDS(counts_files))

        # Move rownames to a column
        counts_reprocessed <- counts_reprocessed %>%
            tibble::rownames_to_column("Sample")

        # Clean the Sample names
        counts_reprocessed$Sample <- counts_reprocessed$Sample %>%
            sub("filtered_", "", .) %>%
            sub("_\\d+\\.fastq", "", .)

        # Move cleaned names back to rownames
        counts_reprocessed <- counts_reprocessed %>%
            tibble::column_to_rownames("Sample")

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
                tax_reprocessed2$BioProject <- study_prefix
                tax_reprocessed$Sequence <- rownames(tax_reprocessed)
            } else {
                stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
            }
            
            } else {
            tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
        }


        # ----- Unzip and read taxonomy -----
        unzipped = unzip(tax_zip, exdir = temp_dir, overwrite = TRUE)
        tax_files = unzipped[grepl("_taxa\\.rds$", unzipped)][1]
        if (length(tax_files) == 0) stop("No *_taxa.rds file found in: ", tax_zip)
        tax_reprocessed <- as.data.frame(readRDS(tax_files))
        tax_reprocessed <- make_taxa_label(tax_reprocessed)

        # Taxonomy rownames = ASVs/Features: prefix if needed
        tax_reprocessed$BioProject <- study_prefix
        tax_reprocessed$Sequence <- rownames(tax_reprocessed)

        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
            counts_reprocessed = aligned$reprocessed
            
            # Skip processing if no rows in counts_reprocessed
            if (nrow(counts_reprocessed) > 0) {
                matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), tax_reprocessed$Sequence)]
                colnames(counts_reprocessed) <- matched_taxa
                counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
                original_names <- colnames(counts_reprocessed)
                counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
            }

            counts_reprocessed2 = aligned$reprocessed
            
            # Skip processing if no rows in counts_reprocessed
            if (nrow(counts_reprocessed2) > 0) {
                matched_taxa <- tax_reprocessed2$Taxa[match(colnames(counts_reprocessed2), tax_reprocessed2$Sequence)]
                colnames(counts_reprocessed2) <- matched_taxa
                counts_reprocessed2 <- collapse_duplicate_columns_exact(counts_reprocessed2)
                original_names2 <- colnames(counts_reprocessed2)
                counts_reprocessed2 <- as.data.frame(lapply(counts_reprocessed2, as.numeric), row.names = rownames(counts_reprocessed2), col.names = original_names2, check.names = FALSE)
                proportions_reprocessed2 <- sweep(counts_reprocessed2, 1, rowSums(counts_reprocessed2), '/')
                all_counts2[[i]] <- counts_reprocessed2
                all_props2[[i]]  <- proportions_reprocessed2
                all_taxa2[[i]]   <- tax_reprocessed2

            }
        }

        if (nrow(counts_reprocessed) > 0) {
            # ----- Proportions -----
            proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
            
            # Store results only if we have valid data
            all_counts[[i]] <- counts_reprocessed
            all_props[[i]]  <- proportions_reprocessed
            all_taxa[[i]]   <- tax_reprocessed
        } else {
            warning("No valid data found in ", counts_zip)
        }

        cleanup_tempfiles(temp_dir)
        }
    }

    # ----- Merge all dataframes -----
    combined_counts <- all_counts %>%
        imap(~ as.data.frame(.x) %>% rownames_to_column("Sample")) %>%
        reduce(full_join, by = "Sample") %>%
        column_to_rownames("Sample") %>%
        replace(is.na(.), 0)
    combined_props  <- all_props %>%
        imap(~ as.data.frame(.x) %>% rownames_to_column("Sample")) %>%
        reduce(full_join, by = "Sample") %>%
        column_to_rownames("Sample") %>%
        replace(is.na(.), 0)
    combined_taxa = bind_rows(all_taxa) 

    if (!file.exists(file.path(local, "rdp16classified.csv.zip"))) {
        combined_counts2 = all_counts2 %>%
        imap(~ as.data.frame(.x) %>% rownames_to_column("Sample")) %>%
        reduce(full_join, by = "Sample") %>%
        column_to_rownames("Sample") %>%
        replace(is.na(.), 0)
        combined_props2 = all_props2 %>%
        imap(~ as.data.frame(.x) %>% rownames_to_column("Sample")) %>%
        reduce(full_join, by = "Sample") %>%
        column_to_rownames("Sample") %>%
        replace(is.na(.), 0)
        combined_taxa2 = bind_rows(all_taxa2)
        write.csv(combined_taxa2, file = file.path(local, "rdp16classified.csv"), row.names = TRUE)
    } else {
        combined_counts2 = all_counts2 %>%
        imap(~ as.data.frame(.x) %>% rownames_to_column("Sample")) %>%
        reduce(full_join, by = "Sample") %>%
        column_to_rownames("Sample") %>%
        replace(is.na(.), 0)
        combined_props2 = all_props2 %>%
        imap(~ as.data.frame(.x) %>% rownames_to_column("Sample")) %>%
        reduce(full_join, by = "Sample") %>%
        column_to_rownames("Sample") %>%
        replace(is.na(.), 0)
        combined_taxa2 = read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
    }

    if (!raw) {
        counts = fill_na_zero_numeric(counts)
        proportions = fill_na_zero_numeric(proportions)
        combined_counts = fill_na_zero_numeric(combined_counts)
        combined_props = fill_na_zero_numeric(combined_props)
        combined_counts2 = fill_na_zero_numeric(combined_counts2)
        combined_props2 = fill_na_zero_numeric(combined_props2)
    }

    # ---- return ----
    return(list(
    scale=scale, 
    metadata=metadata, 
    counts=list(
      original=counts, 
      reprocessed=list(rdp19 = combined_counts, rdp16 = combined_counts2)
    ),
    tax=list(
      original=tax,
      reprocessed=list(rdp19 = combined_taxa, rdp16 = combined_taxa2)
    ),
    proportions=list(
      original=proportions,
      reprocessed=list(rdp19 = combined_props, rdp16 = combined_props2)
    )
  ))
}