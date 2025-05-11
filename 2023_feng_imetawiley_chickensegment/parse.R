parse_2023_feng_imetawiley_chickensegment <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }
    if (!is.logical(raw) || length(raw) != 1) {
        stop("`raw` must be a single logical value (TRUE or FALSE)")
    }
    if (!is.logical(align) || length(align) != 1) {
        stop("`align` must be a single logical value (TRUE or FALSE)")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2023_feng_imetawiley_chickensegment")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJNA817429_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA817429_dada2_taxa.rds.zip")
    counts_16s_zip       <- file.path(local, "16S.csv.zip")
    counts_ITS_zip       <- file.path(local, "ITS.csv.zip")
    metadata_SRA_zip     <- file.path(local, "SraRunTable (34).csv.zip")
    scalereadin          <- file.path(local, "scale.csv.zip")

    # ----- scale and metadata -----
    sra = read_zipped_table(metadata_SRA_zip, row.names=NULL) %>% 
            rename(Accession = Run)
    sra <- sra %>%
        mutate(`Sample Name` = gsub("_", "E", `Sample Name`, fixed = TRUE)) %>%
        separate(`Sample Name`, into = c("Sample", "Amplicontype"), sep = "E(?=[^E]*$)")
    scale = read_zipped_table(scalereadin, row.names=NULL)
    metadata = left_join(sra, scale %>% select(c("Segment","Date", "Number", "Sample")), by = "Sample")
    scale = scale %>% select(-c("Segment","Date", "Number")) %>% 
        mutate(log2_qPCR_16S = ifelse(10^qPCR_log10_16S > 0, log2(10^qPCR_log10_16S), NA)) %>%
        mutate(log2_qPCR_ITS= ifelse(10^qPCR_log10_ITS > 0, log2(10^qPCR_log10_ITS), NA)) %>% 
        rename(log10_qPCR_16S = qPCR_log10_16S, log10_qPCR_ITS = qPCR_log10_ITS)

    # ----- counts, tax, proportions -----
    counts_16s = read_zipped_table(counts_16s_zip, row.names=NULL) %>%
                    select(-c("Group1","Group2")) %>% rename(Sample = index) %>% 
                    as.data.frame()
    rownames(counts_16s) = counts_16s$Sample
    counts_16s$Sample = NULL

    counts_ITS = read_zipped_table(counts_ITS_zip, row.names=NULL) %>%
                    select(-c("Group1","Group2")) %>% rename(Sample = index) %>% 
                    as.data.frame()
    rownames(counts_ITS) = counts_ITS$Sample
    counts_ITS$Sample = NULL

    if (!raw) {
        aligned = rename_and_align(counts_original = counts_16s, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
        counts_16s = aligned$counts_original
        original_names <- colnames(counts_16s)
        counts_16s <- as.data.frame(lapply(counts_16s, as.numeric), row.names = rownames(counts_16s), col.names = original_names, check.names = FALSE)
        aligned = rename_and_align(counts_original = counts_ITS, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
        counts_ITS = aligned$counts_original
        original_names <- colnames(counts_ITS)
        counts_ITS <- as.data.frame(lapply(counts_ITS, as.numeric), row.names = rownames(counts_ITS), col.names = original_names, check.names = FALSE)
    }

    tax_16s <- data.frame(taxonomy = colnames(counts_16s), stringsAsFactors = FALSE)
    tax_ITS <- data.frame(taxonomy = colnames(counts_ITS), stringsAsFactors = FALSE)

    if (any(grepl("^D_\\d+__", tax_16s$taxonomy))) {
    tax_16s <- tax_16s %>%
        mutate(original_taxonomy = taxonomy) %>%
        separate(taxonomy,
                into = c("D0", "D1", "D2", "D3", "D4", "D5", "D6"),
                sep = ";", fill = "right") %>%
        transmute(
        taxonomy = original_taxonomy,
        Kingdom = gsub("^D_0__", "", D0),
        Phylum  = gsub("^D_1__", "", D1),
        Class   = gsub("^D_2__", "", D2),
        Order   = gsub("^D_3__", "", D3),
        Family  = gsub("^D_4__", "", D4),
        Genus   = gsub("^D_5__", "", D5),
        Species = gsub("^D_6__", "", D6)
        )
    }
    tax_16s <- tax_16s %>%
        mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), NA, .)))
    tax_16s[is.na(tax_16s)] <- "unclassified"
    tax_16s = make_taxa_label(tax_16s)
    rownames(tax_16s) <- tax_16s$taxonomy
    if (!raw) {
        matched_taxa <- tax_16s$Taxa[match(colnames(counts_16s), rownames(tax_16s))]
        colnames(counts_16s) <- matched_taxa
        counts_16s <- as.data.frame(t(rowsum(t(counts_16s), group = colnames(counts_16s))))
    }
    tax_ITS$ogtaxonomy <- tax_ITS$taxonomy
    tax_ITS <- tax_ITS %>%
    separate(taxonomy,
            into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
            sep = ";", fill = "right") %>%
    mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))
    tax_ITS <- tax_ITS %>%
        mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), NA, .)))
    tax_ITS[is.na(tax_ITS)] <- "unclassified"
    tax_ITS = make_taxa_label(tax_ITS)
    rownames(tax_ITS) <- tax_ITS$ogtaxonomy
    if (!raw) {
        matched_taxa <- tax_ITS$Taxa[match(colnames(counts_ITS), rownames(tax_ITS))]
        colnames(counts_ITS) <- matched_taxa
        counts_ITS <- collapse_duplicate_columns_exact(counts_ITS)
        original_names <- colnames(counts_ITS)
        counts_ITS <- as.data.frame(lapply(counts_ITS, as.numeric), row.names = rownames(counts_ITS), col.names = original_names, check.names = FALSE)
    }

    # --- Compute proportions from counts ---
    proportions_16s <- sweep(counts_16s, MARGIN = 1,STATS  = rowSums(counts_16s), FUN = "/")
    proportions_ITS <- sweep(counts_ITS, MARGIN = 1,STATS  = rowSums(counts_ITS), FUN = "/")

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
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
        aligned = rename_and_align(counts_original = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
        counts_reprocessed = aligned$counts_original
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed = collapse_duplicate_columns_exact(counts_reprocessed)
        original_names <- colnames(counts_reprocessed)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
    }

    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')

    if (!raw) {
        counts_16s = fill_na_zero_numeric(counts_16s)
        counts_ITS = fill_na_zero_numeric(counts_ITS)
        proportions_16s = fill_na_zero_numeric(proportions_16s)
        proportions_ITS = fill_na_zero_numeric(proportions_ITS)
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = list(
                sixteenS = counts_16s,
                ITS = counts_ITS
            ),
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = list(
                sixteenS = proportions_16s,
                ITS = proportions_ITS
            ),
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = list(
                sixteenS = tax_16s,
                ITS = tax_ITS
            ),
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))


}