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
    temp_dir <- tempfile("metadata")
    dir.create(temp_dir)
    unzip(metadata_zip, exdir = temp_dir)
    xlsx_file <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)[1]
    metadata4 <- list(
        sheet1 = read_excel(xlsx_file, sheet = 1),
        sheet4 = read_excel(xlsx_file, sheet = 4)
    )

    metadata4$sheet1 = metadata4$sheet4 %>% rename(`Library Name` = `BIOSAMPLE`)
    
    # Merge metadata and metadata2 by Library Name
    metadata1 <- full_join(metadata2, metadata4$sheet1, 
                     by = "Library Name")
    
    # Merge with metadata4 sheet1 by fastq ID
    metadata <- full_join(metadata, metadata4$sheet1,
                     by = c("Library Name" = "fastq ID"))
    
    metadata = bind_rows(metadata1, metadata) %>% rename(Sample = `Library Name`)
    
    scale = metadata %>% select(c("Sample", "16S rRNA gene copies / DNA ng (miseq)", "16S rRNA gene copies  / DNA ng (ddPCR)")) %>% 
        mutate(log2_miseq = ifelse(`16S rRNA gene copies / DNA ng (miseq)` > 0, log2(`16S rRNA gene copies / DNA ng (miseq)`), NA)) %>% 
        mutate(log2_ddPCR = ifelse(`16S rRNA gene copies  / DNA ng (ddPCR)` > 0, log2(`16S rRNA gene copies  / DNA ng (ddPCR)`), NA)) %>% 
        mutate(log10_miseq = ifelse(`16S rRNA gene copies / DNA ng (miseq)` > 0, log10(`16S rRNA gene copies / DNA ng (miseq)`), NA)) %>% 
        mutate(log10_ddPCR = ifelse(`16S rRNA gene copies  / DNA ng (ddPCR)` > 0, log10(`16S rRNA gene copies  / DNA ng (ddPCR)`), NA))

    cleanup_tempfiles(temp_dir)

    # ---- Reprocessed data -----
    all_counts <- list()
    all_props  <- list()
    all_taxa   <- list()

    if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {

        for (i in seq_along(repro_counts_rds_zip)) {
            counts_zip <- repro_counts_rds_zip[i]
            tax_zip    <- repro_tax_zip[i]

            # ----- Study prefix -----
            study_prefix <- gsub("_dada2_counts\\.rds\\.zip$", "", basename(counts_zip))

            # ----- Unzip and read counts -----
            temp_dir <- tempfile("repro")
            dir.create(temp_dir)
            unzipped <- unzip(counts_zip, exdir = temp_dir, overwrite = TRUE)
            counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
            if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
            counts_reprocessed <- as.data.frame(readRDS(counts_file))

            # ----- Unzip and read taxonomy -----
            unzipped = unzip(tax_zip, exdir = temp_dir, overwrite = TRUE)
            tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
            if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
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

            cleanup_tempfiles(temp_dir)
        }
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