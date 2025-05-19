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
    all_counts2 <- list()
    all_props2  <- list()
    all_taxa2   <- list()       

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

            # ----- rdp16 -----
            if (!file.exists(file.path(local,"rdp16classified.csv.zip"))) {
            if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
                required_pkgs <- c("dada2", "Biostrings")
                missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
                if (length(missing_pkgs) > 0) {
                    stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                        ". Please install them before running this function.")
                }
                seqs <- Biostrings::DNAStringSet(colnames(counts_reprocessed))
                rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
                tax_reprocessed2 = make_taxa_label(rdpclassified) 
                tax_reprocessed2$BioProject <- study_prefix
                tax_reprocessed2$Sequence <- rownames(tax_reprocessed2)
                } else {
                stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
            }
            
            } else {
                tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
            }

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
                counts_reprocessed2 = aligned$reprocessed
                matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
                matched_taxa2 <- tax_reprocessed2$Taxa[match(colnames(counts_reprocessed2), rownames(tax_reprocessed2))]
                colnames(counts_reprocessed) <- matched_taxa
                colnames(counts_reprocessed2) <- matched_taxa2
                counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
                counts_reprocessed2 <- collapse_duplicate_columns_exact(counts_reprocessed2)
                original_names <- colnames(counts_reprocessed)
                original_names2 <- colnames(counts_reprocessed2)
                counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
                counts_reprocessed2 <- as.data.frame(lapply(counts_reprocessed2, as.numeric), row.names = rownames(counts_reprocessed2), col.names = original_names2, check.names = FALSE)
                proportions_reprocessed2 <- sweep(counts_reprocessed2, 1, rowSums(counts_reprocessed2), '/')
            }

            # ----- Proportions -----
            proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')

            # Store results
            all_counts[[i]] <- counts_reprocessed
            all_props[[i]]  <- proportions_reprocessed
            all_taxa[[i]]   <- tax_reprocessed

            all_counts2[[i]] <- counts_reprocessed2
            all_props2[[i]]  <- proportions_reprocessed2
            if (!file.exists(file.path(local, "rdp16classified.csv.zip"))) {
                all_taxa2[[i]]   <- tax_reprocessed2
            }   

            cleanup_tempfiles(temp_dir)
        }
    }

    # ----- Merge all dataframes -----
    combined_counts <- bind_rows(all_counts)
    combined_props  <- bind_rows(all_props)
    combined_taxa   <- bind_rows(all_taxa)

    if (!file.exists(file.path(local, "rdp16classified.csv.zip"))) {
        combined_counts2 <- bind_rows(all_counts2)
        combined_props2  <- bind_rows(all_props2)
        combined_taxa2   <- bind_rows(all_taxa2)
    } else {
        combined_counts2 <- NA
        combined_props2  <- NA
        combined_taxa2   <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
    }

    if (!raw) {
        counts = fill_na_zero_numeric(counts)
        proportions = fill_na_zero_numeric(proportions)
        combined_counts = fill_na_zero_numeric(combined_counts)
        combined_props = fill_na_zero_numeric(combined_props)
        combined_counts2 = fill_na_zero_numeric(combined_counts2)
        combined_props2 = fill_na_zero_numeric(combined_props2)
    }

    return(list(
        scale=scale, 
        metadata=metadata, 
        counts=list(
            original=counts, 
            reprocessed=list(rdp19 = combined_counts, rdp16 = combined_counts2)
        ),
        tax=list(
            original=tax,
            reprocessed=list(rdp19 = combined_taxa, rdp16 = combined_taxa2)
        ),
        proportions=list(
            original=proportions,
            reprocessed=list(rdp19 = combined_props, rdp16 = combined_props2)
        )
    ))
}