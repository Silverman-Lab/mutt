parse_2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct <- function() {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    local               <- "2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct/"
    repro_counts_zips   <- c(
                            paste0(local, "PRJNA394877_dada2_merged_nochim.rds.zip"),
                            paste0(local, "PRJNA548153_dada2_taxonomy_merged.rds.zip"),
                            paste0(local, "PRJNA606262_dada2_taxonomy_merged.rds.zip"),
                            paste0(local, "PRJNA607574_dada2_taxonomy_merged.rds.zip"),
                            paste0(local, "PRJNA545312_dada2_merged_nochim.rds.zip")
    )
    repro_tax_zips      <- c(
                            paste0(local, "PRJNA394877_dada2_taxonomy_merged.rds.zip"),
                            paste0(local, "PRJNA548153_dada2_taxonomy_merged.rds.zip"),
                            paste0(local, "PRJNA606262_dada2_taxonomy_merged.rds.zip"),
                            paste0(local, "PRJNA607574_dada2_taxonomy_merged.rds.zip"),
                            paste0(local, "PRJNA545312_dada2_taxonomy_merged.rds.zip")
    )
    scale_zip           <- paste0(local, "Liao2021_scale.csv")
    metadata_zip        <- paste0(local, "Liao_2021_metadata.csv")
    counts_zip          <- paste0(local, "Liao_2021_16S.csv")

    read_zipped_csv <- function(zip_path) {
        if (file.exists(zip_path)) {
            csv_file <- unzip(zip_path, list = TRUE)$Name[1]
            read.csv(unz(zip_path, csv_file), row.names = 1, check.names = FALSE)
        } else {
            warning(paste("File not found:", zip_path))
            return(NA)
        }
    }

    # ------ original counts ------
    counts_original <- read_zipped_csv(counts_16s_zip)

    if (!is.na(counts_original)[1]) {
        original_taxa <- colnames(counts_original)
        taxon_ids <- paste0("Taxon_", seq_len(nrow(counts_original)))
        colnames(counts_original) <- taxon_ids

        # Create taxa mapping data frame
        tax_original_mice <- data.frame(
        Taxon = taxon_ids,
        Original_Taxa = original_taxa,
        stringsAsFactors = FALSE
        )

        # ------ proportions from counts ------
        proportions_original <- t(counts_original_mice)
        proportions_original[] <- lapply(
        proportions_original,
        function(col) {
            if (is.numeric(col)) {
            total <- sum(col, na.rm = TRUE)
            if (total == 0) return(rep(NA, length(col)))
            return(col / total)
            } else {
            return(col)
            }
        }
        )
    } else {
        proportions_original <- NA
        tax_original_mice <- NA
    }

    # ---- scale and metadata -----
    scale     <- read_zipped_csv(scale_16s_zip)
    metadata  <- read_zipped_csv(metadata_16s_zip)

    for (i in seq_along(repro_counts_zips)) {
        # Counts
        temp_rds <- tempfile(fileext = ".rds")
        unzip(repro_counts_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
        rds_file <- list.files(dirname(temp_rds), pattern = "\\.rds$", full.names = TRUE)[1]
        seqtab_nochim <- readRDS(rds_file)
        rpt_mat <- t(seqtab_nochim)
        counts <- as.data.frame(rpt_mat)
        counts$Sequence <- rownames(counts)
        counts <- counts[, c("Sequence", setdiff(names(counts), "Sequence"))]
        rownames(counts) <- paste0("Taxon_", seq_len(nrow(counts)))
        counts_reprocessed_list[[i]] <- counts

        # Taxonomy
        temp_tax <- tempfile(fileext = ".rds")
        unzip(repro_tax_zips[i], exdir = dirname(temp_tax), overwrite = TRUE)
        tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
        taxonomy_matrix <- readRDS(tax_file)
        rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
        tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
        tax_reprocessed_list[[i]] <- tax_table
    }

    # Combine reprocessed
    counts_reprocessed <- bind_rows(counts_reprocessed_list)
    proportions_reprocessed <- counts_reprocessed
    proportions_reprocessed[-1] <- lapply(
        counts_reprocessed[-1],
        function(col) col / sum(col)
    )
    tax_reprocessed <- bind_rows(tax_reprocessed_list)

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