parse_2023_galla_scientificreports_spikeinmockvariousorganismssampletypes <- function(raw = FALSE, align = FALSE) {
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
    local <- file.path("2023_galla_scientificreports_spikeinmockvariousorganismssampletypes")

    # ----- File paths -----
    repro_counts_rds_zip <- c(
        file.path(local, "PRJNA703791_dada2_counts.rds.zip"),
        file.path(local, "PRJNA734187_dada2_counts.rds.zip")
    )
    repro_tax_zip        <- c(
        file.path(local, "PRJNA703791_dada2_taxa.rds.zip"),
        file.path(local, "PRJNA734187_dada2_taxa.rds.zip")
    )

    metadata_SRA1_zip    <- file.path(local, "SraRunTable (38).csv.zip")
    metadata_SRA2_zip    <- file.path(local, "SraRunTable (39).csv.zip")
    metadata_zip         <- file.path(local, "41598_2023_30916_MOESM3_ESM.xlsx.zip")
    metadata1_zip        <- file.path(local, "table1.csv.zip")

    # ----- counts, tax, proportions -----
    counts = NA
    tax = NA
    proportions = NA

    # ----- scale and metadata -----
    metadata = read_zipped_table(metadata_SRA1_zip, row.names = NULL) %>% rename(Accession = Run)
    metadata2 = read_zipped_table(metadata_SRA2_zip, row.names = NULL) %>% rename(Accession = Run)
    metadata3 = read_zipped_table(metadata1_zip, row.names = NULL) 
    temp_dir <- tempdir()
    unzip(metadata_zip, exdir = temp_dir)
    xlsx_file <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)[1]
    metadata4 <- list(
        sheet1 = read_excel(xlsx_file, sheet = 1),
        sheet4 = read_excel(xlsx_file, sheet = 4)
    )

    metadata4$sheet1 = metadata4$sheet4 %>% rename(`Library Name` = `BIOSAMPLE`)
    
    # Merge metadata and metadata2 by Library Name
    metadata1 <- merge(metadata2, metadata4$sheet1, 
                     by = "Library Name",
                     all = TRUE)
    
    # Merge with metadata4 sheet1 by fastq ID
    metadata <- merge(metadata, metadata4$sheet1,
                     by.x = "Library Name",
                     by.y = "fastq ID", 
                     all = TRUE)
    
    metadata = bind_rows(metadata1, metadata) %>% rename(Sample = `Library Name`)
    
    scale = metadata %>% select(c("Sample", "16S rRNA gene copies / DNA ng (miseq)", "16S rRNA gene copies  / DNA ng (ddPCR)")) %>% 
        mutate(log2_miseq = ifelse(`16S rRNA gene copies / DNA ng (miseq)` > 0, log2(`16S rRNA gene copies / DNA ng (miseq)`), NA)) %>% 
        mutate(log2_ddPCR = ifelse(`16S rRNA gene copies  / DNA ng (ddPCR)` > 0, log2(`16S rRNA gene copies  / DNA ng (ddPCR)`), NA)) %>% 
        mutate(log10_miseq = ifelse(`16S rRNA gene copies / DNA ng (miseq)` > 0, log10(`16S rRNA gene copies / DNA ng (miseq)`), NA)) %>% 
        mutate(log10_ddPCR = ifelse(`16S rRNA gene copies  / DNA ng (ddPCR)` > 0, log10(`16S rRNA gene copies  / DNA ng (ddPCR)`), NA))

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

        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
            counts_reprocessed = aligned$reprocessed
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), tax_reprocessed$Sequence)]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
            original_names <- colnames(counts_reprocessed)
            counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
        }

        # ----- Proportions -----
        proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')

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