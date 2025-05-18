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
    scale <- read_zipped_table(scale_zip, row.names = NULL) %>%
            mutate(
                log2_qPCR = qPCR16S * log2(10) 
            ) %>%
            rename(log10_qPCR = qPCR16S)
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
            temp_rds <- tempfile("repro")
            dir.create(temp_rds)
            unzipped = unzip(repro_counts_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
            counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
            if (is.na(counts_file)) stop(paste("No *_counts.rds file found for index", i))
            counts_reprocessed <- as.data.frame(readRDS(counts_file))

            # Unzip and read taxonomy
            unzipped = unzip(repro_tax_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
            tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
            if (is.na(tax_file)) stop(paste("No *_taxa.rds file found for index", i))
            tax_reprocessed <- as.data.frame(readRDS(tax_file))
            tax_reprocessed <- make_taxa_label(tax_reprocessed)

            # Replace sequence columns with taxon labels
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))

            if (!raw) {
                aligned <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
                counts_reprocessed <- aligned$counts_reprocessed
                original_names <- colnames(counts_reprocessed)
                counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
            }
            # proportions
            proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')

            # Label with study name based on zip filename prefix
            study_id <- sub("_.*$", "", basename(tools::file_path_sans_ext(repro_counts_zips[i])))
            counts_reprocessed$Study <- study_id
            proportions_reprocessed$Study <- study_id
            tax_reprocessed$Study <- study_id

            counts_reprocessed_list[[i]] <- counts_reprocessed
            proportions_reprocessed_list[[i]] <- proportions_reprocessed
            tax_reprocessed_list[[i]] <- tax_reprocessed

            cleanup_tempfiles(temp_rds)
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