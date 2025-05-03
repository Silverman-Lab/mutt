parse_2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct <- function(raw = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
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


    read_zipped_table <- function(zip_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE) {
      if (file.exists(zip_path)) {
        inner_file <- unzip(zip_path, list = TRUE)$Name[1]
        con <- unz(zip_path, inner_file)
        read.table(con, sep = sep, header = header, row.names = row.names, check.names = check.names, stringsAsFactors = FALSE)
      } else {
        warning(paste("File not found:", zip_path))
        return(NA)
      }
    }

    # ------ original counts ------
    counts_original <- read_zipped_table(counts_16s_zip)

    if (!is.na(counts_original)[1]) {
        original_taxa <- colnames(counts_original)

        # Create taxa mapping data frame
        tax_original_mice <- data.frame(
        Taxa = original_taxa,
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
    scale     <- read_zipped_table(scale_16s_zip)
    metadata  <- read_zipped_table(metadata_16s_zip)

    # Helper to create lowest-rank taxon label
    make_taxa_label <- function(df) {
        tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
        prefixes  <- c("k", "p", "c", "o", "f", "g")
        if (!all(tax_ranks %in% colnames(df))) {
            stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
        }
        df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
            x[is.na(x) | trimws(x) == ""] <- "unclassified"
            return(x)
        })
        df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
            for (i in length(tax_ranks):1) {
                if (tax_row[i] != "unclassified") {
                    return(paste0(prefixes[i], "_", tax_row[i]))
                }
            }
            return("unclassified")
        })
        return(df)
    }

    # Process multiple zipped RDS files
    counts_reprocessed_list <- list()
    proportions_reprocessed_list <- list()
    tax_reprocessed_list <- list()

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

        # Compute proportions
        proportions <- counts
        proportions[] <- lapply(counts, function(col) col / sum(col))

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