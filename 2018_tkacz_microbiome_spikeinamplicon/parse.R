parse_2018_tkacz_microbiome_spikeinamplicon <- function(raw = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
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












    read_zipped_table <- function(zip_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE) {
    if (file.exists(zip_path)) {
        inner_file <- unzip(zip_path, list = TRUE)$Name[1]
        con <- unz(zip_path, inner_file)
        read.table(con, sep = sep, header = header, row.names = row.names, check.names = check.names, stringsAsFactors = FALSE)
    } else {
        warning(paste("File not found:", zip_path))
        return(NULL)
    }
    }

    # ---- Taxa labeling function ----
    make_taxa_label <- function(df) {
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    prefixes <- c("k", "p", "c", "o", "f", "g")
    df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
        x[is.na(x) | trimws(x) == ""] <- "unclassified"
        x
    })
    df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
        if (tax_row["Genus"] != "unclassified") {
        return(paste0("g_", tax_row["Genus"]))
        }
        for (i in (length(tax_ranks) - 1):1) {
        if (tax_row[i] != "unclassified") {
            return(paste0("uc_", prefixes[i], "_", tax_row[i]))
        }
        }
        return("unclassified")
    })
    return(df)
    }

    # ---- Process and merge all studies ----
    all_counts <- list()
    all_tax <- list()
    all_metadata <- list()

    for (i in seq_along(repro_counts_rds_zip)) {
    # Unzip and read RDS
    temp_dir <- tempfile("study")
    dir.create(temp_dir)
    unzip(repro_counts_rds_zip[i], exdir = temp_dir, overwrite = TRUE)
    counts_file <- list.files(temp_dir, pattern = "_counts\\.rds$", full.names = TRUE)
    stopifnot(length(counts_file) == 1)
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
    counts <- as.data.frame(t(rowsum(t(counts), group = colnames(counts))))

    # Store
    all_counts[[i]] <- counts
    all_tax[[i]] <- tax

    # Read SRA metadata
    metadata <- read_zipped_table(metadata_SRA_zip[i], sep = ",", header = TRUE, row.names = 1)
    all_metadata[[i]] <- metadata
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

    # ---- Proportions ----
    proportions_merged <- counts_merged
    proportions_merged[] <- lapply(counts_merged, function(col) col / sum(col))


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