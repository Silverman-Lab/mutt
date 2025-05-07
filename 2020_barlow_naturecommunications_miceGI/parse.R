parse_2020_barlow_naturecommunications_miceGI <- function(raw = FALSE) {
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
    local <- file.path("2020_barlow_naturecommunications_miceGI")

    # ----- File paths -----
    sra_metadata_zip     <- file.path(local, "SraRunTable (32).csv.zip")
    counts               <- NA  # No original counts provided
    sraplusload_zip      <- file.path(local, "SraRunTable_total_loads.csv.zip")
    proportions_zip      <- file.path(local, "Relative_Abundance_Table.csv.zip")
    scale_zip            <- file.path(local, "absolute_abundance_calculated.csv.zip")
    repro_counts_rds_zip <- file.path(local, "PRJNA575097_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA575097_dada2_taxa.rds.zip")

    # ---- helper functions ----
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
        if (missing(x)) return(NULL)
        if (is.data.frame(x)) {
            x[] <- lapply(x, function(y) if (is.numeric(y)) replace(y, is.na(y), 0) else y)
        } else if (is.matrix(x) && is.numeric(x)) {
            x[is.na(x)] <- 0
        } else if (is.list(x)) {
            x <- lapply(x, fill_na_zero_numeric)
        }
        x
    }
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

    # --- proportions ---
    prop_csv <- unzip(proportions_zip, list = TRUE)$Name[1]
    prop_path <- unzip(proportions_zip, files = prop_csv, exdir = tempdir(), overwrite = TRUE)
    proportions <- read.csv(prop_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

    # metadatabackup <- proportions %>%
    #     rownames_to_column("Sample") %>%
    #     select(c(Sample, Diet, Site, Day, mouse, Cage)) %>% 
    #     mutate(Sample = str_replace_all(Sample, "_", " "))

    proportions <- proportions %>% select(-c(Diet, Site, Day, mouse, Cage))

    # --- Assign Taxon_1 ... Taxon_N as rownames and store original names ---
    original_taxa <- colnames(proportions)

    # --- Create taxa dataframe (lookup table) ---
    tax <- data.frame(taxonomy = original_taxa, stringsAsFactors = FALSE)
    rownames(tax) <- tax$taxonomy
    if (any(grepl("^D_\\d+__", tax$taxonomy))) {
    tax <- tax %>%
        separate(taxonomy,
                into = paste0("D", 0:6),
                sep = ";", fill = "right") %>%
        mutate(across(everything(), ~ gsub("^D_\\d+__", "", .))) %>%
        mutate(across(everything(), ~ ifelse(. == "" | . == "__", "unclassified", .))) %>%
        transmute(
        Kingdom = D0,
        Phylum  = D1,
        Class   = D2,
        Order   = D3,
        Family  = D4,
        Genus   = D5,
        Species = D6
        )
    }

    tax_original = make_taxa_label(tax) %>% as.data.frame(stringsAsFactors = FALSE)
    rownames(tax_original) <- tax_original$taxonomy
    proportions <- as.data.frame(proportions, stringsAsFactors = FALSE)
    if (!raw) {
        matched_taxa <- tax_original$Taxa[match(colnames(proportions), rownames(tax_original))]
        colnames(proportions) <- matched_taxa
        proportions <- as.data.frame(t(rowsum(t(proportions), group = colnames(proportions))))
    }

    # ----- metadata ------
    sra_metadata = read_zipped_table(sra_metadata_zip, row.names = NULL) %>% 
                    mutate(Sample = str_replace_all(`Sample Name`, "_", " ")) %>% 
                    rename(Accession = Run)
    
    # --- Scale (ddPCR) ---
    scale_csv <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_csv, exdir = tempdir(), overwrite = TRUE)
    scale <- read.csv(scale_path, row.names = NULL, stringsAsFactors = FALSE)
    metadata <- full_join(sra_metadata, scale, by = "Sample") 
    metadata = read_zipped_table(sraplusload_zip, row.names = NULL) %>% rename(Accession = Run) %>%
                    full_join(metadata, by = "Accession")
    dups <- duplicated(as.list(df))
    metadata <- metadata[, !dups] 
    scale = metadata %>% select(c("Sample" ,"Accession", "Total Load (16S Copies/g)"))
    metadata = metadata %>% select(-c("Total Load (16S Copies/g)"))
    
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
      counts_original = fill_na_zero_numeric(counts_original)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_original = fill_na_zero_numeric(proportions_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }

    # ---- return ----
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
        scale = genomicscale
    ))
}
