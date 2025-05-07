parse_2023_kallastu_research_foodscience_food <- function(raw = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }
    # Load needed libraries
    library(tidyverse)
    library(readxl)
    library(readr)

    # -------- local path -----------
    local                       <- file.path("2023_kallastu_research_foodscience_food")

    # -------- file paths -----------
    metadata_zip                <- file.path(local, "SraRunTable.csv.zip")
    repro_counts_rds_zip        <- file.path(local, "PRJNA861123_dada2_counts.rds.zip")
    repro_tax_zip               <- file.path(local, "PRJNA861123_dada2_taxa.rds.zip")

    # ---- helper functions ----

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

    # ------- metadata -------------
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- as.data.frame(read.csv(metadata_path, row.names = NULL, stringsAsFactors = FALSE)) %>%
                    rename(Accession = Run, Sample = Sample.Name)

    # ------- scale ----------------
    ## Scale is directly presented as a table in the paper
    scale_data <- data.frame(
        "20St" = c(1.55e9, NA, 1.56e9, 1.16e9),
        `K6 PMA` = c(1.28e8, 1.53e7, 5.96e7, 1.71e8),
        `K6 TOT` = c(1.86e9, NA, 1.26e8, 3.06e8),
        `K5 PMA` = c(3.82e8, 3.64e7, 7.04e7, 3.73e8),
        `K5 TOT` = c(5.90e8, NA, 1.34e8, 5.66e8),
        `K4 PMA` = c(3.97e8, 1.17e8, 2.61e7, 4.12e8),
        `K4 TOT` = c(8.56e8, NA, 1.49e8, 4.60e8),
        `20St 2.5PMA` = c(1.51e9, NA, NA, NA),
        `20St 1PMA`   = c(1.39e9, NA, NA, NA),
        `20St 0.5PMA` = c(1.74e9, NA, NA, NA)
    )
    rownames(scale_data) <- c("Spike_In", "Plating_cfu_per_g", "FC", "qPCR")
    colnames(scale_data) <- c(
        "20St", "K6 PMA", "K6 TOT", "K5 PMA", "K5 TOT",
        "K4 PMA", "K4 TOT", "20St 2.5PMA", "20St 1PMA", "20St 0.5PMA"
    )
    scale = as.data.frame(t(scale_data)) %>% rownames_to_column("Sample")

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

    if (!raw) {
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }
    
    return(list(
        scale=scale, 
        metadata=metadata, 
        counts=list(
            original = NA, 
            reprocessed = counts_reprocessed
        ),
        proportions=list(
            original = NA,
            reprocessed = proportions_reprocessed
        ),
        tax=list(
            original = NA,
            reprocessed = tax_reprocessed
        )
        )
    )
}