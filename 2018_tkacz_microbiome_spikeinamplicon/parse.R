parse_2018_tkacz_microbiome_spikeinamplicon <- function(raw = FALSE, align = FALSE) {
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
    local <- file.path("2018_tkacz_microbiome_spikeinamplicon")

    # ----- File paths -----
    repro_counts_rds_zip <- c(
        file.path(local, "PRJEB22042_dada2_counts.rds.zip"),
        file.path(local, "PRJEB22043_dada2_counts.rds.zip")
    )
    repro_tax_zip        <- c(
        file.path(local, "PRJEB22042_dada2_taxa.rds.zip"),
        file.path(local, "PRJEB22043_dada2_taxa.rds.zip")
    )
    metadata_SRA_zip     <- c(
        file.path(local, "SraRunTable (38).csv.zip"),
        file.path(local, "SraRunTable (39).csv.zip")
    )

    # --- initialize ---
    counts = NA
    tax = NA
    proportions = NA

    # ---- Process and merge all studies ----
    all_counts <- list()
    all_tax <- list()
    all_metadata <- list()

    if (all(file.exists(repro_counts_rds_zip))) {
        for (i in seq_along(repro_counts_rds_zip)) {
        # Unzip and read RDS
        temp_dir <- tempfile("study")
        dir.create(temp_dir)
        unzipped = unzip(repro_counts_rds_zip[i], exdir = temp_dir, overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
        counts <- as.data.frame(readRDS(counts_file))

        unzip(repro_tax_zip[i], exdir = temp_dir, overwrite = TRUE)
        tax_file <- list.files(temp_dir, pattern = "_taxa\\.rds$", full.names = TRUE)
        stopifnot(length(tax_file) == 1)
        tax <- as.data.frame(readRDS(tax_file))

        # Label taxa
        tax <- make_taxa_label(tax)

        # Match column names to tax labels
        matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
        colnames(counts) <- matched_taxa
        counts <- collapse_duplicate_columns_exact(counts)
        original_names <- colnames(counts)
        counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)

        # Store
        all_counts[[i]] <- counts
        all_tax[[i]] <- tax

        # Read SRA metadata
        metadata <- read_zipped_table(metadata_SRA_zip[i], sep = ",", header = TRUE, row.names = 1)
        all_metadata[[i]] <- metadata
        cleanup_tempfiles(temp_dir)
        }
    }

    # ---- Merge ----
    counts_merged <- Reduce(function(a, b) {
    full_join(a, b, by = "row.names") %>%
        mutate(Row.names = coalesce(!!sym("Row.names"), !!sym("Taxa"))) %>%
        column_to_rownames("Row.names") %>%
        mutate_all(~ replace_na(., 0))
    }, all_counts)

    tax_merged <- bind_rows(all_tax)
    metadata_merged <- bind_rows(all_metadata)

    if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_merged, metadata = metadata_merged, scale = scale, by_col = "sampleID", align = align, study_name = basename(local))
        counts_merged <- aligned$reprocessed
        original_names <- colnames(counts_merged)
        counts_merged <- as.data.frame(lapply(counts_merged, as.numeric), row.names = rownames(counts_merged), col.names = original_names, check.names = FALSE)
    }

    # ---- Proportions ----
    proportions_merged <- sweep(counts_merged, 1, rowSums(counts_merged), '/')

    if (!raw) {
      counts = fill_na_zero_numeric(counts)
      counts_merged = fill_na_zero_numeric(counts_merged)
      proportions = fill_na_zero_numeric(proportions)
      proportions_merged = fill_na_zero_numeric(proportions_merged)
    }

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts,
            reprocessed = counts_merged
        ),
        proportions = list(
            original = proportions,
            reprocessed = proportions_merged
        ),
        tax = list(
            original = tax,
            reprocessed = tax_merged
        ),
        scale = scale,
        metadata = metadata_merged
    ))
}