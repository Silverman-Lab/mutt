parse_2023_galla_scientificreports_spikeinmockvariousorganismssampletypes <- function(raw=FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2023_galla_scientificreports_spikeinmockvariousorganismssampletypes")

    # ----- File paths -----
    repro_counts_rds_zip <- c(
        file.path(local, "PRJNA703791_dada2_counts.rds.zip")
        file.path(local, "PRJNA734187_dada2_counts.rds.zip")
    )
    repro_tax_zip        <- c(
        file.path(local, "PRJNA703791_dada2_taxa.rds.zip")
        file.path(local, "PRJNA734187_dada2_taxa.rds.zip")
    )

    metadata_SRA1_zip    <- file.path(local, "SraRunTable (38).csv.zip")
    metadata_SRA2_zip    <- file.path(local, "SraRunTable (39).csv.zip")
    metadata_zip         <- file.path(local, "41598_2023_30916_MOESM3_ESM.xlsx.zip")
    metadata1_zip        <- file.path(local, "table1.csv.zip")

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
    fill_na_zero_numeric <- function(x) {
        if (is.data.frame(x)) {
            x[] <- lapply(x, function(y) if (is.numeric(y)) replace(y, is.na(y), 0) else y)
        } else if (is.matrix(x) && is.numeric(x)) {
            x[is.na(x)] <- 0
        } else if (is.list(x)) {
            x <- lapply(x, fill_na_zero_numeric)
        }
        x
    }

    # ----- counts, tax, proportions -----

    # ----- scale and metadata -----


    # ---- Reprocessed data -----
    all_counts <- list()
    all_props  <- list()
    all_taxa   <- list()

    for (i in seq_along(repro_counts_zips)) {
        counts_zip <- repro_counts_zips[i]
        tax_zip    <- repro_tax_zips[i]

        # ----- Study prefix -----
        study_prefix <- gsub("_dada2_counts\\.rds\\.zip$", "", basename(counts_zip))

        # ----- Unzip and read counts -----
        temp_rds <- tempfile(fileext = ".rds")
        unzip(counts_zip, exdir = dirname(temp_rds), overwrite = TRUE)
        rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
        if (length(rds_files) == 0) stop("No *_counts.rds file found in: ", counts_zip)
        counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

        # ----- Unzip and read taxonomy -----
        temp_tax <- tempfile(fileext = ".rds")
        unzip(tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
        tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
        if (length(tax_files) == 0) stop("No *_taxa.rds file found in: ", tax_zip)
        tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))
        tax_reprocessed <- make_taxa_label(tax_reprocessed)

        # Taxonomy rownames = ASVs/Features: prefix if needed
        tax_reprocessed$BioProject <- study_prefix
        tax_reprocessed$Sequence <- rownames(tax_reprocessed)

        # Rename counts columns using taxonomy
        if (!raw) {
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), tax_reprocessed$Sequence)]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
        }

        # ----- Proportions -----
        proportions_reprocessed <- counts_reprocessed
        proportions_reprocessed[] <- lapply(
        proportions_reprocessed,
        function(col) col / sum(col)
        )

        # Store results
        all_counts[[i]] <- counts_reprocessed
        all_props[[i]]  <- proportions_reprocessed
        all_taxa[[i]]   <- tax_reprocessed
    }

    # ----- Merge all dataframes -----
    combined_counts <- do.call(rbind, all_counts)
    combined_props  <- do.call(rbind, all_props)
    combined_taxa   <- bind_rows(all_taxa)

    if (!raw) {
        counts = fill_na_zero_numeric(counts)
        proportions = fill_na_zero_numeric(proportions)
        combined_counts = fill_na_zero_numeric(combined_counts)
        combined_props = fill_na_zero_numeric(combined_props)
    }

    return(list(
        scale=scale, 
        metadata=metadata, 
        counts=list(
            original=counts, 
            reprocessed=combined_counts
        ),
        tax=list(
            original=tax,
            reprocessed=combined_taxa
        ),
        proportions=list(
            original=proportions,
            reprocessed=combined_props
        )
    ))
}