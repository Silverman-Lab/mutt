parse_2025_wagner_frontiersinmiccrobiology_flowpiglets <- function(raw = FALSE, align = FALSE) {
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
    local <- file.path("2025_wagner_frontiersinmiccrobiology_flowpiglets")

    # ----- File paths -----
    sra_metadata_zip     <- file.path(local, "SraRunTable (23).csv.zip")
    sra_metadata_zip_2   <- file.path(local, "SraRunTable (24).csv.zip")
    counts               <- file.path(local, "counts.csv.zip")
    scale_zip            <- file.path(local, "scale.csv.zip")
    repro_counts_rds_zip <- c(
        file.path(local, "PRJNA1229264_dada2_counts.rds.zip"),
        file.path(local, "PRJNA800240_dada2_counts.rds.zip")
    )
    repro_tax_zip        <- c(
        file.path(local, "PRJNA1229264_dada2_taxa.rds.zip"),
        file.path(local, "PRJNA800240_dada2_taxa.rds.zip")
    )

    convert_scientific <- function(x) {
    # Extract base and exponent from format like "9.37 × 107"
    matches <- regexec("([0-9.]+)\\s*×\\s*10([0-9]+)", x)
    parts <- regmatches(x, matches)
    
    sapply(parts, function(p) {
        if (length(p) == 3) {
        base <- as.numeric(p[2])
        exp  <- as.numeric(p[3])
        base * 10^exp
        } else {
        NA_real_
        }
    })
    }

    # ---- scale ----
    scale <- read_zipped_table(scale_zip, row.names = NULL) %>%
    as.data.frame() %>%
    mutate(Sample = paste0("swine", `Swine no.`, "_d", day)) %>%
    relocate(Sample, .before = 1) %>%
    select(-`Swine no.`, -day) %>%
    mutate(
        cells_per_gram_feces_mean = convert_scientific(cells_per_gram_feces_mean),
        cells_per_gram_feces_sd   = convert_scientific(cells_per_gram_feces_sd)
    ) %>% 
    mutate(log2_FC_cells_g_mean  = ifelse(cells_per_gram_feces_mean > 0, log2(cells_per_gram_feces_mean), NA),
            log10_FC_cells_g_mean = ifelse(cells_per_gram_feces_mean > 0, log10(cells_per_gram_feces_mean), NA),
            log2_FC_cells_g_sd    = ifelse(cells_per_gram_feces_sd   > 0, log2(cells_per_gram_feces_sd), NA),
            log10_FC_cells_g_sd   = ifelse(cells_per_gram_feces_sd   > 0, log10(cells_per_gram_feces_sd), NA))


    # ---- metadata ----
    metadata = read_zipped_table(sra_metadata_zip, row.names = NULL) %>% as.data.frame() %>% rename(Accession = Run) %>%
                mutate(Sample = gsub("^swine_(\\d+)_day_(\\d+)$", "swine\\1_d\\2", host_subject_id, perl = TRUE)) %>% relocate(Sample, .before = 1)


    # ---- counts, tax ----
    counts_original = read_zipped_table(counts) %>% t() %>% as.data.frame(stringsAsFactors = FALSE)
    counts_original[] <- lapply(counts_original, function(x) as.numeric(as.character(x)))
    tax_original = data.frame(taxonomy = colnames(counts_original))

    # ---- proportions ----
    proportions_original <- sweep(counts_original, MARGIN = 1,STATS  = rowSums(counts_original), FUN = "/")

    if (!raw) {
        aligned = rename_and_align(counts_original = counts_original,  proportions_original = proportions_original,
                                        metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
            counts_original = aligned$counts_original
            proportions_original = aligned$proportions_original 
            original_names <- colnames(counts_original)
            counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
            original_names <- colnames(proportions_original)
            proportions_original <- as.data.frame(lapply(proportions_original, as.numeric), row.names = rownames(proportions_original), col.names = original_names, check.names = FALSE)
    }

    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
    tax_reprocessed2 <- NA
    proportions_reprocessed2 <- NA
    counts_reprocessed2 <- NA

    if (all(file.exists(repro_counts_rds_zip))) {
        # ----- Reprocessed counts from RDS ZIP -----
        temp_dir <- tempfile("repro")
        dir.create(temp_dir)
        unzipped <- unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
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
        unzipped <- unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))

        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa ----- 
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed,
                                        metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
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
        counts_original = fill_na_zero_numeric(counts_original)
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_original = fill_na_zero_numeric(proportions_original)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
        counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
    }

    # ---- return ----
    return(list(
        counts = list(
            original = counts_original,
            reprocessed = list(counts = counts_reprocessed, proportions = proportions_reprocessed)
        ),
        tax = list(
            original = tax_original,
            reprocessed = list(tax = tax_reprocessed, proportions = proportions_reprocessed2)
        ),
        proportions = list(
            original = proportions_original,
            reprocessed = list(proportions = proportions_reprocessed, proportions2 = proportions_reprocessed2)
        ),
        metadata = metadata,
        scale = scale
    ))
}