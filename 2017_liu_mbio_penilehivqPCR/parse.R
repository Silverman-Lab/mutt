parse_2017_liu_mbio_penilehivqPCR <- function(raw = FALSE, align = FALSE) {
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
    local <- file.path("2017_liu_mbio_penilehivqPCR")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJNA1233249_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA1233249_dada2_taxa.rds.zip")
    scale_16s_zip        <- file.path(local, "PMID-28743816_samples-v1.csv.zip")
    metadata_16s_zip     <- file.path(local, "metadata.csv.zip")

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA

    # ---- scale and metadata -----
    scale     <- read_zipped_table(scale_16s_zip, row.names = NULL)
    metadata  <- read_zipped_table(metadata_16s_zip, row.names= NULL) 

    scale <- scale %>%
        select(ID, `16S Qty Mean`) %>%
        rename(Sample_name = ID) %>%
        left_join(metadata %>% select(Sample_name, Accession), by = "Sample_name") %>%
        mutate(log2_qPCR = ifelse(`16S Qty Mean`>0, log2(`16S Qty Mean`), NA)) %>%
        mutate(log10_qPCR = ifelse(`16S Qty Mean`>0, log10(`16S Qty Mean`), NA))

    # ------ original counts ------

    # READ IN ORIGINAL COUNTS -------------------------- Needs update

    # ---- tax and proportions -----
    if (!is.na(counts_original)[1]) {

        original_taxa <- colnames(counts_original)

        # Create taxa mapping data frame
        tax_original <- data.frame(
        Taxa = original_taxa,
        stringsAsFactors = FALSE
        )

        # ------ proportions from counts ------
        proportions_original <- sweep(counts_original, 1, rowSums(counts_original, na.rm = TRUE), FUN = "/")

    } else {
        proportions_original <- NA
        tax_original <- NA
    }

    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
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
            align <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample_name", align = align, study_name = basename(local))
            counts_reprocessed = align$reprocessed
        }

        # taxa
        if (!raw) {
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
        }

        # proportions reprocessedF_cell_ml
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