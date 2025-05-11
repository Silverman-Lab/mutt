parse_2020_zemb_microOpen_spike <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }
    if (!is.logical(raw) || length(raw) != 1) {
        stop("`raw` must be a single logical value (TRUE or FALSE)")
    }
    if (!is.logical(align) || length(align) != 1) {
        stop("`align` must be a single logical value (TRUE or FALSE)")
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
    sra_metadata <- read.csv(sra_metadata_path, row.names = NULL, stringsAsFactors = FALSE) %>% rename(Accession = Run)
    sra_metadata$Sample.Name <- paste0("oz1802-", str_extract(sra_metadata$Sample.Name, "(?<=Tube).*")) 
    rownames(sra_metadata) <- sra_metadata$Sample.Name

    # --- Scale (qPCR) ---
    scale <- read_zipped_table(scale_zip, row.names = 1) %>% rownames_to_column("Sample.Name") %>%
        rename(
            standard_measure1_copies_per_ul = "Std_M1(compies/ulPCR)",
            standard_measure2_copies_per_ul = "Std_M2(compies/ulPCR)", 
            standard_measure3_copies_per_ul = "Std_M3(compies/ulPCR)",
            standard_copies_in_tube = "Std_in_tube (copies)",
            pcr_efficiency = "efficiency",
            total_16s_measure1_copies_per_ul = "total16S_M1(copies/ul pcr)",
            total_16s_measure2_copies_per_ul = "total16S_M2(copies/ul pcr)", 
            total_16s_measure3_copies_per_ul = "total16S_M3(copies/ul pcr)",
            average_16s_copies_per_ul = "average16S(copies/ul pcr)",
            average_16s_copies_in_tube = "average16S(copies_in_tube)",
            average_16s_copies_per_mg = "average16S(copies/mg)"
        ) %>%
        mutate(
            log2_copies_ul = ifelse(average_16s_copies_per_ul > 0, log2(average_16s_copies_per_ul), NA),
            log10_copies_ul = ifelse(average_16s_copies_per_ul > 0, log10(average_16s_copies_per_ul), NA),
            log2_copies_mg = ifelse(average_16s_copies_per_mg > 0, log2(average_16s_copies_per_mg), NA),
            log10_copies_mg = ifelse(average_16s_copies_per_mg > 0, log10(average_16s_copies_per_mg), NA),
            log2_copies_tube = ifelse(average_16s_copies_in_tube > 0, log2(average_16s_copies_in_tube), NA),
            log10_copies_tube = ifelse(average_16s_copies_in_tube > 0, log10(average_16s_copies_in_tube), NA)
        )


    # Merge sra_metadata and metadata
    metadata = cbind(metadata, sra_metadata[rownames(metadata), , drop = FALSE]) %>% rename(SampleID = Sample.Name) %>% rownames_to_column("Sample.Name")


    # --- Counts ---
    counts_csv <- unzip(counts_zip, list = TRUE)$Name[1]
    counts_path <- unzip(counts_zip, files = counts_csv, exdir = tempdir(), overwrite = TRUE)
    counts = read.csv(counts_path, row.names = 1, stringsAsFactors = FALSE)
    metadata2 <- counts %>% rownames_to_column("Sample.Name") %>% select(c("Sample.Name", "X.1"))
    counts <- counts %>% select(-c("X.1"))
    colnames(counts) <- sapply(colnames(counts), function(x) {
      tube_num <- str_match(x, "Tube(\\d+)")[,2]
      if (!is.na(tube_num)) {
        paste0("oz1802-", tube_num)
      } else {
        x
      }
    })
    counts <- counts %>% t() %>% as.data.frame()

    # --- Create taxa dataframe ---
    original_taxa <- colnames(counts)
    tax <- data.frame(
    Taxa = original_taxa,
    stringsAsFactors = FALSE
    )

    if (!raw) {
        aligned <- rename_and_align(counts_original = counts, metadata = metadata, scale = scale, by_col = "Sample.Name", align = align, study_name = basename(local))
        counts <- aligned$counts_original
        original_names <- colnames(counts)
        counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
    }

    # --- Compute proportions from counts ---
    proportions <- sweep(counts, 1, rowSums(counts, na.rm = TRUE), '/')

    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(repro_counts_rds_zip))) {
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
        if (!raw) {
            aligned <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample.Name", align = align, study_name = basename(local))
            counts_reprocessed <- aligned$reprocessed
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
            original_names <- colnames(counts_reprocessed)
            counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
        }

        # proportions reprocessed
        proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
    }

    if (!raw) {
      counts = fill_na_zero_numeric(counts)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions = fill_na_zero_numeric(proportions)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }

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
