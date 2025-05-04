parse_2022_cvandevelde_ismecommunications_culturedflowhumanfecal <- function(raw = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
            stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2022_cvandevelde_ismecommunications_culturedflowhumanfecal")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJEB51873_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJEB51873_dada2_taxa.rds.zip")
    scale_16s_zip        <- file.path(local, "VandeVelde2022_scale.csv.zip")
    counts_16s_zip       <- file.path(local, "VandeVelde_2022_16S.csv.zip")
    metadata_16s_zip     <- file.path(local, "VandeVelde_2022_metadata.csv.zip")

    read_zipped_csv <- function(zip_path) {
        if (file.exists(zip_path)) {
            csv_file <- unzip(zip_path, list = TRUE)$Name[1]
            read.csv(unz(zip_path, csv_file), row.names = 1, check.names = FALSE)
        } else {
            warning(paste("File not found:", zip_path))
            return(NA)
        }
    }

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA

    # ------ original counts ------
    counts_original <- read_zipped_csv(counts_16s_zip)

    if (!is.na(counts_original)[1]) {
        original_taxa <- colnames(counts_original)
        taxon_ids <- paste0("Taxon_", seq_len(nrow(counts_original)))
        colnames(counts_original) <- taxon_ids

        # Create taxa mapping data frame
        tax_original <- data.frame(
        Taxon = taxon_ids,
        Original_Taxa = original_taxa,
        stringsAsFactors = FALSE
        )

        # ------ proportions from counts ------
        proportions_original <- t(counts_original)
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
        tax_original <- NA
    }

    # ---- scale and metadata -----
    scale     <- read_zipped_csv(scale_16s_zip)
    metadata  <- read_zipped_csv(metadata_16s_zip)


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
    make_taxa_label <- function(df) {
        tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
        prefixes  <- c("k", "p", "c", "o", "f", "g")
        if (!all(tax_ranks %in% colnames(df))) {
            stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
        }
        df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
            x[is.na(x) | trimws(x) == ""] <- "unclassified"
            x
        })
        df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
            if (tax_row["Genus"] != "unclassified") {
            return(paste0("g_", tax_row["Genus"]))
            }
            for (i in (length(tax_ranks)-1):1) {  # skip Genus
            if (tax_row[i] != "unclassified") {
                return(paste0("uc_", prefixes[i], "_", tax_row[i]))
            }
            }
            return("unclassified")
        })
        return(df)
    }
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    # accessions to sampleIDs is study specific: IF NEED BE

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
