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
    counts_reprocessed2 <- NA
    proportions_reprocessed2 <- NA
    tax_reprocessed <- NA
    tax_reprocessed2 <- NA
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
    if (file.exists(repro_counts_rds_zip)) {
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
        tax_reprocessed = make_taxa_label(tax_reprocessed)
        
        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
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

        if (!raw) {
            counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
            proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
            counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
            proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
        }
    } else {
        counts_reprocessed = NA
        proportions_reprocessed = NA
        tax_reprocessed = NA
        counts_reprocessed2 = NA
        proportions_reprocessed2 = NA
        tax_reprocessed2 = NA
    }

    if (!raw) {
        counts_original = fill_na_zero_numeric(counts_original)
        proportions_original = fill_na_zero_numeric(proportions_original)
    }

    cleanup_tempfiles(temp_dir)

    metadata <- metadata %>%
      # Convert empty strings to NA
      mutate(across(everything(), ~ifelse(. == "" | . == " ", NA, .))) %>%
      # Convert specific columns to factors
      mutate(across(c("Prefix", "Well"), factor)) %>%
      # Clean up timepoint values and convert to numeric
      mutate(Timepoint = str_replace_all(Timepoint, "h", ""),
             Timepoint = as.numeric(Timepoint)) %>%
      # Convert numeric columns
      mutate(across(where(is.character), ~ifelse(grepl("^[0-9.]+$", .), as.numeric(.), .)))


    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original,
            reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
        ),
        proportions = list(
            original = proportions_original,
            reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
        ),
        tax = list(
            original = tax_original,
            reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
        ),
        scale = scale,
        metadata = metadata
    ))
}
