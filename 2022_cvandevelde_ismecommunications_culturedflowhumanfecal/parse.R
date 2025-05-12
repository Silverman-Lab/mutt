parse_2022_cvandevelde_ismecommunications_culturedflowhumanfecal <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
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

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2022_cvandevelde_ismecommunications_culturedflowhumanfecal")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJEB51873_dada2_counts.rds.zip") 
    repro_tax_zip        <- file.path(local, "PRJEB51873_dada2_taxa.rds.zip") 
    scale_zip            <- file.path(local, "VandeVelde2022_scale.csv.zip")
    counts_zip           <- file.path(local, "VandeVelde_2022_16S.csv.zip")
    metadata_zip         <- file.path(local, "VandeVelde_2022_metadata.csv.zip")
    sra_zip              <- file.path(local, "SraRunTable (39).csv.zip")

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
    scale <- NA
    metadata <- NA
    sra <- NA

    # ---- scale and metadata -----
    tryCatch({
        scale <- read_zipped_table(scale_zip, row.names = NULL) %>% 
            as.data.frame() 
        scale <- scale %>% 
            rename(Sample = sampleID) %>%
            mutate(log2_FC = ifelse(10^load > 0, log2(10^load), NA)) %>%
            rename(log10_FC = load)
    }, error = function(e) {
        warning("Error reading scale file: ", e$message)
    })

    # no way to connect the metadata to the scale and counts because non similar sampleIDs, I suspect it might have to do with plate names but cant figure out.
    # need to email authors about this, but also, there are a ton more samples in the metadata and reprocessed data and in the flow cytometry data that needs to be added.
    # Actually the counts might be in supplementary data, but I dont really know if theyre useable for this -- i dont think so. so need to review study and ask authors about the genomic
    # data connection to our scale. 
    tryCatch({
        metadata <- read_zipped_csv(metadata_zip) %>% 
            as.data.frame() 
        metadata <- metadata %>% 
            rename(Sample = Sample_ID)
    }, error = function(e) {
        warning("Error reading metadata file: ", e$message)
    })

    tryCatch({
        sra <- read_zipped_table(sra_zip) %>% 
            rename(Accession = Run)
    }, error = function(e) {
        warning("Error reading SRA file: ", e$message)
    })

    # ------ original counts ------
    tryCatch({
        counts_original <- read_zipped_csv(counts_zip)

        if (!is.na(counts_original)[1]) {
            original_taxa <- colnames(counts_original)

            # Create taxa mapping data frame
            tax_original <- data.frame(
                Taxa = original_taxa,
                stringsAsFactors = FALSE
            )

            if (!raw) {
                aligned = rename_and_align(counts_original = counts_original, 
                                        metadata = metadata, 
                                        scale = scale, 
                                        by_col = "Sample", 
                                        align = align, 
                                        study_name = basename(local))
                counts_original = aligned$counts_original
                original_names <- colnames(counts_original)
                counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
            }

            # ------ proportions from counts ------
            proportions_original <- sweep(counts_original, 1, rowSums(counts_original), '/')
        }
    }, error = function(e) {
        warning("Error processing original counts: ", e$message)
    })

    # ----- Reprocessed counts from RDS ZIP -----
    tryCatch({
        temp_rds <- tempfile("repro")
        dir.create(temp_rds)
        unzipped = unzip(repro_counts_rds_zip, exdir = temp_rds, overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(counts_file))

        # ----- Taxonomy reprocessed -----
        unzipped = unzip(repro_tax_zip, exdir = temp_rds, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))
        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)
    }, error = function(e) {
        warning("Error processing reprocessed taxonomy: ", e$message)
    })

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw && !is.na(counts_reprocessed)[1]) {
        tryCatch({
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed, 
                                    metadata = metadata, 
                                    scale = scale, 
                                    by_col = "Sample_name", 
                                    align = align, 
                                    study_name = basename(local))
            counts_reprocessed = aligned$reprocessed
            asv_seqs <- trimws(colnames(counts_reprocessed))
            tax_seqs <- trimws(rownames(tax_reprocessed))
            seq_to_taxa <- setNames(tax_reprocessed$Taxa, tax_seqs)
            matched_taxa <- seq_to_taxa[asv_seqs]
            matched_taxa[is.na(matched_taxa)] <- "unclassified"
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
            original_names <- colnames(counts_reprocessed)
            counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
        }, error = function(e) {
            warning("Error aligning reprocessed data: ", e$message)
        })
    }

    # proportions reprocessed
    if (!is.na(counts_reprocessed)[1]) {
        tryCatch({
            proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
        }, error = function(e) {
            warning("Error calculating reprocessed proportions: ", e$message)
        })
    }

    if (!raw) {
        tryCatch({
            if (!is.na(counts_original)[1]) {
                counts_original = fill_na_zero_numeric(counts_original)
            }
            if (!is.na(proportions_original)[1]) {
                proportions_original = fill_na_zero_numeric(proportions_original)
            }
            if (!is.na(counts_reprocessed)[1]) {
                counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
            }
            if (!is.na(proportions_reprocessed)[1]) {
                proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
            }
        }, error = function(e) {
            warning("Error filling NA values: ", e$message)
        })
    }

    cleanup_tempfiles(temp_rds)

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original,
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = proportions_original,
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = tax_original,
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))
}
