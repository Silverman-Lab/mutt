parse_2020_zemb_microOpen_spike <- function(raw = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }

    # Load libraries
    library(tidyverse)
    library(readxl)
    library(stringr)
    library(readr)

    # ----- Local base directory -----
    local <- file.path("2020_zemb_microOpen_spike")

    # ----- File paths -----
    counts_zip           <- file.path(local, "zemb_counts.csv.zip")
    metadata_zip         <- file.path(local, "zemb_metadata.csv.zip")
    sra_metadata_zip     <- file.path(local, "SraRunTable.csv.zip")
    scale_zip            <- file.path(local, "zemb_qPCR.csv.zip")
    repro_counts_rds_zip <- file.path(local, "PRJNA531076_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA531076_dada2_taxa.rds.zip")

    # --- Metadata ---
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- read.csv(metadata_path, row.names = 1, stringsAsFactors = FALSE)

    sra_metadata_csv <- unzip(sra_metadata_zip, list = TRUE)$Name[1]
    sra_metadata_path <- unzip(sra_metadata_zip, files = sra_metadata_csv, exdir = tempdir(), overwrite = TRUE)
    sra_metadata <- read.csv(sra_metadata_path, row.names = 1, stringsAsFactors = FALSE)
    rownames(sra_metadata) <- sra_metadata$Sample.Name
    sra_metadata <- subset(sra_metadata, select = -Sample.Name)
    rownames(sra_metadata) <- paste0("oz1802-", str_extract(rownames(sra_metadata), "(?<=Tube).*")) 

    # --- Counts ---
    counts_csv <- unzip(counts_zip, list = TRUE)$Name[1]
    counts_path <- unzip(counts_zip, files = counts_csv, exdir = tempdir(), overwrite = TRUE)
    counts = read.csv(counts_path, row.names = 1, stringsAsFactors = FALSE)
    colnames(counts) <- paste0("oz1802-", str_extract(colnames(counts), "(?<=Tube).*"))

    # --- Scale (qPCR) ---
    scale_csv <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_csv, exdir = tempdir(), overwrite = TRUE)
    scale <- read.csv(scale_path, row.names = 1, stringsAsFactors = FALSE)

    # counts <- counts[, colnames(counts) %in% rownames(metadata)]
    # metadata <- metadata[rownames(metadata) %in% colnames(counts), ]
    # scale = scale[rownames(scale) %in% colnames(counts), ]
    # sra_metadata <- sra_metadata[rownames(sra_metadata) %in% colnames(counts), ]

    # Merge sra_metadata and metadata
    metadata = cbind(metadata, sra_metadata[rownames(metadata), , drop = FALSE])

    # Check alignment
    # all.equal(rownames(metadata), rownames(scale))
    # all.equal(colnames(counts), rownames(scale))

    # --- Create taxa dataframe ---
    original_taxa <- rownames(counts)
    tax <- data.frame(
    Taxa = original_taxa,
    stringsAsFactors = FALSE
    )

    # --- Compute proportions from counts ---
    proportions <- counts
    proportions[] <- lapply(
    proportions,
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

    # Return structured list
    return(list(
        counts = list(
            original = counts,
            reprocessed = counts_reprocessed
        ),
        tax = list(
            original = tax,
            reprocessed = tax_reprocessed
        ),
        proportions = list(
            original = proportions,
            reprocessed = proportions_reprocessed
        ),
        metadata = metadata,
        scale = scale
    ))
}
