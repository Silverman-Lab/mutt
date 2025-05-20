parse_2023_kallastu_research_foodscience_food <- function(raw = FALSE, align = FALSE) {
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
    local                       <- file.path("2023_kallastu_research_foodscience_food")

    # -------- file paths -----------
    metadata_zip                <- file.path(local, "SraRunTable.csv.zip")
    repro_counts_rds_zip        <- file.path(local, "PRJNA861123_dada2_counts.rds.zip")
    repro_tax_zip               <- file.path(local, "PRJNA861123_dada2_taxa.rds.zip")

    # ------- metadata -------------
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- as.data.frame(read.csv(metadata_path, row.names = NULL, stringsAsFactors = FALSE)) %>%
                    rename(Accession = Run, Sample = Sample.Name)

    # ------- scale ----------------
    ## Scale is directly presented as a table in the paper
    scale_data <- data.frame(
        "20St" = c(1.55e9, NA, 1.56e9, 1.16e9),
        `K6 PMA` = c(1.28e8, 1.53e7, 5.96e7, 1.71e8),
        `K6 TOT` = c(1.86e9, NA, 1.26e8, 3.06e8),
        `K5 PMA` = c(3.82e8, 3.64e7, 7.04e7, 3.73e8),
        `K5 TOT` = c(5.90e8, NA, 1.34e8, 5.66e8),
        `K4 PMA` = c(3.97e8, 1.17e8, 2.61e7, 4.12e8),
        `K4 TOT` = c(8.56e8, NA, 1.49e8, 4.60e8),
        `20St 2.5PMA` = c(1.51e9, NA, NA, NA),
        `20St 1PMA`   = c(1.39e9, NA, NA, NA),
        `20St 0.5PMA` = c(1.74e9, NA, NA, NA)
    )
    rownames(scale_data) <- c("Spike_In", "Plating_cfu_per_g", "FC", "qPCR")
    colnames(scale_data) <- c(
        "20St", "K6 PMA", "K6 TOT", "K5 PMA", "K5 TOT",
        "K4 PMA", "K4 TOT", "20St 2.5PMA", "20St 1PMA", "20St 0.5PMA"
    )
    scale = as.data.frame(t(scale_data)) %>% rownames_to_column("Sample") %>% 
                    mutate(log2_FC = ifelse(FC > 0, log2(FC), NA)) %>% 
                    mutate(log10_FC = ifelse(FC > 0, log10(FC), NA)) %>% 
                    mutate(log2_Plating_cfu_per_g = ifelse(Plating_cfu_per_g > 0, log2(Plating_cfu_per_g), NA)) %>% 
                    mutate(log10_Plating_cfu_per_g = ifelse(Plating_cfu_per_g > 0, log10(Plating_cfu_per_g), NA)) %>% 
                    mutate(log2_qPCR = ifelse(qPCR > 0, log2(qPCR), NA)) %>% 
                    mutate(log10_qPCR = ifelse(qPCR > 0, log10(qPCR), NA))

    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
        temp_dir <- tempfile("repro")
        dir.create(temp_dir)
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

        # ----- Taxonomy reprocessed ----
        unzipped = unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))

        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
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
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }
    
    return(list(
        scale=scale, 
        metadata=metadata, 
        counts=list(
            original = NA, 
            reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed)
        ),
        proportions=list(
            original = NA,
            reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed)
        ),
        tax=list(
            original = NA,
            reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed)
        )
        )
    )
}