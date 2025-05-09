parse_2019_ji_nature_human_soil <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
            stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                ". Please install them before running this function.")
    }
    if (!is.logical(align)) {
      stop("align must be a logical value")
    }
    if (!is.logical(raw)) {
      stop("raw must be a logical value")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2019_ji_nature_human_soil")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJNA541083_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA541083_dada2_taxa.rds.zip")
    metadata_zip         <- file.path(local, "metadata/Gut_Metadata.csv.zip")
    sra_zip              <- file.path(local, "metadata/SraRunTable (38).csv.zip")

    # ---- original counts, tax, and proportions ----
    counts_original = NA
    tax_original = NA
    proportions_original = NA
    

    # ---- scale andmetadata ----
    sra = read_zipped_table(sra_zip, row.names= NULL) %>% rename(Accession = Run, sampleID = `sample name`)
    metadata = read_zipped_table(metadata_zip, row.names= NULL) 
    metadata = merge(metadata, sra, by = "sampleID")
    scale = metadata %>% select(sampleID, `gDNA quant ng/uL`, `PCR quant`) %>%
                        mutate(log2_gDNA_quant = ifelse(`gDNA quant ng/uL` > 0, log2(`gDNA quant ng/uL`), NA)) %>%
                        mutate(log10_gDNA_quant = ifelse(`gDNA quant ng/uL` > 0, log10(`gDNA quant ng/uL`), NA)) %>%
                        mutate(log2_PCR_quant = ifelse(`PCR quant` > 0, log2(`PCR quant`), NA)) %>%
                        mutate(log10_PCR_quant = ifelse(`PCR quant` > 0, log10(`PCR quant`), NA))

    if (all(file.exists(repro_counts_rds_zip))) {
        # ----- Reprocessed counts from RDS ZIP -----
        temp_rds <- tempfile(fileext = ".rds")
        unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

        rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
        if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

        # ----- Taxonomy reprocessed -----
        temp_tax <- tempfile(fileext = ".rds")
        unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

        tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
        if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            align = rename_and_align(counts_reprocessed = counts_reprocessed, metadata, scale, by_col = "sampleID", align = align, study_name = basename(local))
            counts_reprocessed <- align$reprocessed
        }

        # taxa
        if (!raw) {
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
        }

        # proportions reprocessed
        proportions_reprocessed = counts_reprocessed
        proportions_reprocessed[-1] <- lapply(
            counts_reprocessed[-1],
            function(col) col / sum(col)
        )
    }

    if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_original = fill_na_zero_numeric(proportions_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }

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