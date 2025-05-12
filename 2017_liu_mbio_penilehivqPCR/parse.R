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
    sra_zip              <- file.path(local, "SraRunTable (40).csv.zip")

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA

    # ---- scale and metadata -----
    scale        <- read_zipped_table(scale_16s_zip, row.names = NULL)
    metadata16s  <- read_zipped_table(metadata_16s_zip, row.names= NULL) %>% rename(Sample = Accession, Sample_ID = `Sample ID`)
    sra          <- read_zipped_table(sra_zip, row.names = NULL) %>% rename(Accession = Run)

    metadata <- metadata16s %>%
        left_join(sra, by = "Sample_ID") 

    scale <- scale %>%
        select(ID, `16S Qty Mean`) %>%
        rename(Sample_name = ID) %>%
        left_join(metadata %>% select(Sample_name), by = "Sample_name") %>%
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

        # convert to numeric
        original_names <- colnames(counts_original)
        counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)

        # ------ proportions from counts ------
        proportions_original <- sweep(counts_original, 1, rowSums(counts_original, na.rm = TRUE), FUN = "/")

    } else {
        counts_original <- NA
        proportions_original <- NA
        tax_original <- NA
    }

    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
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

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            aligned <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample_name", align = align, study_name = basename(local))
            counts_reprocessed = aligned$reprocessed
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
            original_names <- colnames(counts_reprocessed)
            counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
        }

        # proportions reprocessed
        proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
        cleanup_tempfiles(temp_rds)
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