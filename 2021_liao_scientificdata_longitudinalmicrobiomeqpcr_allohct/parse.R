parse_2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct <- function(raw = FALSE, align = FALSE) {
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
    local <- file.path("2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct")

    # ----- File paths -----
    repro_counts_zips <- c(
    file.path(local, "PRJNA394877_dada2_counts.rds.zip"),
    file.path(local, "PRJNA548153_dada2_counts.rds.zip"),
    file.path(local, "PRJNA606262_dada2_counts.rds.zip"),
    file.path(local, "PRJNA607574_dada2_counts.rds.zip"),
    file.path(local, "PRJNA545312_dada2_counts.rds.zip")
    )

    repro_tax_zips <- c(
    file.path(local, "PRJNA394877_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA548153_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA606262_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA607574_dada2_taxa.rds.zip"),
    file.path(local, "PRJNA545312_dada2_taxa.rds.zip")
    )

    scale_zip    <- file.path(local, "Liao2021_scale.csv.zip")
    metadata_zip <- file.path(local, "Liao_2021_metadata.csv.zip")
    counts_zip   <- file.path(local, "Liao_2021_16S.csv.zip")

    # ---- scale and metadata -----
    scale     <- read_zipped_table(scale_zip, row.names = NULL) %>% mutate(log2_qPCR = ifelse(10^qPCR16S > 0, log2(10^qPCR16S), NA)) %>% rename(log10_qPCR = qPCR16S)
    metadata  <- read_zipped_table(metadata_zip, row.names = NULL)

    # ------ original counts ------
    counts_original <- read_zipped_table(counts_zip)

    if (!is.na(counts_original)[1]) {
        original_taxa <- colnames(counts_original)

        # Create taxa mapping data frame
        tax_original <- data.frame(
        Taxa = original_taxa,
        stringsAsFactors = FALSE
        )

        if (!raw) {
            aligned <- rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
            counts_original <- aligned$counts_original
            original_names <- colnames(counts_original)
            counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
        }
        # proportions
        proportions_original <- sweep(counts_original, MARGIN = 1, STATS = rowSums(counts_original), FUN = "/")
    } else {
        proportions_original <- NA
        tax_original <- NA
    }

    # Process multiple zipped RDS files
    counts_reprocessed_list <- list()
    proportions_reprocessed_list <- list()
    tax_reprocessed_list <- list()

    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA

    if (all(file.exists(repro_counts_zips))) {
        for (i in seq_along(repro_counts_zips)) {
            # Unzip and read counts
            temp_rds <- tempfile(fileext = ".rds")
            unzip(repro_counts_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
            rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
            if (length(rds_files) == 0) stop(paste("No *_counts.rds file found for index", i))
            counts <- as.data.frame(readRDS(rds_files[1]))

            # Unzip and read taxonomy
            temp_tax <- tempfile(fileext = ".rds")
            unzip(repro_tax_zips[i], exdir = dirname(temp_tax), overwrite = TRUE)
            tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
            if (length(tax_files) == 0) stop(paste("No *_taxa.rds file found for index", i))
            tax <- as.data.frame(readRDS(tax_files[1]))
            tax <- make_taxa_label(tax)

            # Replace sequence columns with taxon labels
            matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
            colnames(counts) <- matched_taxa
            counts <- as.data.frame(t(rowsum(t(counts), group = colnames(counts))))

            if (!raw) {
                aligned <- rename_and_align(counts_reprocessed = counts, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
                counts <- aligned$counts_reprocessed
                original_names <- colnames(counts)
                counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
            }
            # proportions
            proportions <- sweep(counts, 1, rowSums(counts), '/')

            # Label with study name based on zip filename prefix
            study_id <- sub("_.*$", "", basename(tools::file_path_sans_ext(repro_counts_zips[i])))
            counts$Study <- study_id
            proportions$Study <- study_id
            tax$Study <- study_id

            counts_reprocessed_list[[i]] <- counts
            proportions_reprocessed_list[[i]] <- proportions
            tax_reprocessed_list[[i]] <- tax
        }

        # Combine all
        counts_reprocessed <- bind_rows(counts_reprocessed_list)
        proportions_reprocessed <- bind_rows(proportions_reprocessed_list)
        tax_reprocessed <- bind_rows(tax_reprocessed_list)
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