parse_2023_pereira_nature_nervous <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "stringr", "Biostrings", "dada2")
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
    library(Biostrings)
    library(dada2)

    # ---------- MANAN PROCESSED BELOW ----------

    # localPath <- file.path("/2023_pereira_nature_nervous/")

    # decompressed_file <- gunzip(paste0(localPath, "all-relevant-data.xlsx.gz"), remove = FALSE)
    
    # totalCounts <- readxl::read_xlsx(decompressed_file, sheet = "Sup. Table 1", skip = 1)
    # totalCounts <- totalCounts[, c("time (hours)", "condition", "medium", "Cells/mL")]
    # totalCounts$index <- paste(totalCounts[["time (hours)"]], totalCounts$condition, totalCounts$medium, sep = "_")
    # totalCounts$replicate <- ave(totalCounts$index, totalCounts$index, FUN = seq_along)
    # totalCounts$index <- paste(totalCounts$index, totalCounts$replicate, sep = "_")
    # totalCounts <- totalCounts[, c("index", "Cells/mL")]
    
    # sampleMetadata <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 2", skip = 3)
    
    # # Step 1: Process sampleMetadata
    # sampleMetadata <- sampleMetadata %>%
    #   mutate(index = paste(`time (hours)`, condition, medium, sep = "_")) %>%
    #   group_by(index) %>%
    #   mutate(replicate = seq_along(index)) %>%
    #   ungroup() %>%
    #   mutate(index = paste(index, replicate, sep = "_"))
    
    # # Step 2: Intersect with totalCounts
    # scale_data <- sampleMetadata %>%
    #   inner_join(totalCounts, by = "index") %>%
    #   dplyr::select("Sequencing sample ID", `Cells/mL`)
    
    # # Step 3: Create scale matrix
    # scale <- scale_data %>%
    #   pivot_wider(names_from = "Sequencing sample ID", values_from = `Cells/mL`) %>%
    #   as.matrix()
    
    # relativeAbundances <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 3", skip = 3)
    # absoluteAbundances <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 4", skip = 3)
    
    # absoluteAbundances <- absoluteAbundances %>%
    #   rename(ASV_ID = `...1`)
    
    # relativeAbundances <- relativeAbundances %>%
    #   rename(ASV_ID = `...1`)
    
    # countsDataFile <- unzip(paste0(localPath, "DADA2_counts_as_matrix_Drugs.txt.zip"))
    
    # countsData <- read.table(countsDataFile, header = TRUE, sep = "\t")
    
    # counts <- countsData %>% column_to_rownames("X") %>% as.matrix()
    
    # proportions <- relativeAbundances %>%
    #   dplyr::select(-c(Phylum, Class, Order, Family, Genus)) %>%
    #   column_to_rownames("ASV_ID") %>%
    #   as.matrix()
    
    # metadata <- sampleMetadata %>%
    #   column_to_rownames("Sequencing sample ID")
    
    # tax <- absoluteAbundances %>%
    #   dplyr::select(ASV_ID, Phylum, Class, Order, Family, Genus) %>%
    #   column_to_rownames("ASV_ID")
    
    # unlink(decompressed_file)
    # unlink(countsDataFile)

    # ------ MAXWELL PROCESSED BELOW ------------

    # ----- Local base directory -----
    local <- file.path("2023_pereira_nature_nervous")

    # ----- File paths -----
    repro_counts_rds_zip                  <- file.path(local, "PRJNA1033532_dada2_counts.rds.zip")
    repro_tax_zip                         <- file.path(local, "PRJNA1033532_dada2_taxa.rds.zip")
    scale_16s_zip                         <- file.path(local, "Pereira2023_scale.csv.zip")
    metadata_16s_zip                      <- file.path(local, "Pereira_2023_metadata.csv.zip")
    DADA2_counts_as_matrix_Drugs_zip      <- file.path(local, "DADA2_counts_as_matrix_Drugs.txt.zip")
    DADA2_counts_as_matrix_FeRecovery_zip <- file.path(local, "DADA2_counts_as_matrix_FeRecovery.txt.zip")
    cleaned_DADA2_ASVs_Drugs_zip          <- file.path(local, "cleaned_DADA2_ASVs_Drugs.fna.zip")
    cleaned_DADA2_ASVs_FeRecovery_zip     <- file.path(local, "cleaned_DADA2_ASVs_Drugs_FeRecovery.fna.zip")
    sra_zip                               <- file.path(local, "SraRunTable (38).csv.zip")

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_original_16 <- NA
    counts_original_19 <- NA
    proportions_original_16 <- NA
    proportions_original_19 <- NA
    tax_original16 <- NA
    tax_original19 <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    counts_reprocessed2 <- NA
    proportions_reprocessed2 <- NA
    tax_reprocessed <- NA
    tax_reprocessed2 <- NA
    scale <- NA
    metadata <- NA

    # ---- scale and metadata -----
    scale         <- read_zipped_table(scale_16s_zip, row.names = NULL) %>% rename(Sample = `Sequencing sample ID`) 
    metadata      <- read_zipped_table(metadata_16s_zip, row.names = NULL) %>% rename(Sample = `Sequencing sample ID`)
    sra          <- read_zipped_table(sra_zip, row.names = NULL) %>% rename(Accession = Run, Sample = `Library Name`) %>% mutate(Sample = gsub("-", "_", Sample))

    metadata      <- full_join(scale, metadata, by = "Sample")
    metadata      <- full_join(sra, metadata, by = "Sample")
    scale         <- scale %>% 
                    dplyr::select(c("Sample", "Cells/mL")) %>% 
                    mutate(`Cells/mL` = as.numeric(`Cells/mL`)) %>%
                    mutate(log2_FC_cells_ml = ifelse(10^`Cells/mL` > 0, log2(10^`Cells/mL`),NA)) %>% 
                    rename(log10_FC_cells_ml = `Cells/mL`)

    # ------ original counts ------
    counts_original1 <- read_zipped_table(DADA2_counts_as_matrix_Drugs_zip, sep = "\t")
    counts_original2 <- read_zipped_table(DADA2_counts_as_matrix_FeRecovery_zip, sep = "\t")

    df1 <- as.data.frame(counts_original1)
    df1$Feature <- rownames(df1)
    df2 <- as.data.frame(counts_original2)
    df2$Feature <- rownames(df2)
    all <- merge(df1, df2, by="Feature", all=TRUE)
    all[is.na(all)] <- 0
    counts_original <- aggregate(. ~ Feature, data=all, FUN=sum)
    rownames(counts_original) <- counts_original$Feature
    counts_original$Feature <- NULL
    counts_original <- as.data.frame(t(counts_original))

    # ------ original tax ------
    fna_file <- unzip(cleaned_DADA2_ASVs_Drugs_zip, list=TRUE)$Name
    tmp <- tempdir()
    unzip(cleaned_DADA2_ASVs_Drugs_zip, files=fna_file, exdir=tmp)
    seqs <- readDNAStringSet(file.path(tmp, fna_file), format="fasta")

    fna_file <- unzip(cleaned_DADA2_ASVs_FeRecovery_zip, list=TRUE)$Name
    tmp <- tempdir()
    unzip(cleaned_DADA2_ASVs_FeRecovery_zip, files=fna_file, exdir=tmp)
    seqs2 <- readDNAStringSet(file.path(tmp, fna_file), format="fasta")

    all_seqs <- c(seqs, seqs2)
    unique_seqs <- unique(all_seqs)

    if (!is.na(counts_original) && !is.na(unique_seqs)) {

        # ----- rdp16 -----    
        if (!file.exists(file.path(local,"rdp16classified_ORIGINAL.csv.zip"))) {
            if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
                rdp16classified <- dada2::assignTaxonomy(unique_seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
                tax_original16 = make_taxa_label(rdp16classified) 
                rownames(tax_original16) <- unique_seqs$names
                write.csv(tax_original16, file = file.path(local, "rdp16classified_ORIGINAL.csv"), row.names = TRUE)
            } else {
                stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
            }
            if (file.exists(file.path("helperdata/rdp_19_toGenus_trainset.fa.gz"))) {
                rdp19classified <- dada2::assignTaxonomy(unique_seqs, file.path("helperdata/rdp_19_toGenus_trainset.fa.gz"), multithread=TRUE) %>% as.data.frame()
                tax_original19 = make_taxa_label(rdp19classified) 
                rownames(tax_original19) <- unique_seqs$names
                write.csv(tax_original19, file = file.path(local, "rdp19classified_ORIGINAL.csv"), row.names = TRUE)
            } else {
                stop("RDP 19 file not detected. please install the helperdata/rdp_19_toGenus_trainset.fa.gz file")
            }
        } else {
            tax_original16 <- read_zipped_table(file.path(local, "rdp16classified_ORIGINAL.csv.zip"), row.names = 1)
            tax_original19 <- read_zipped_table(file.path(local, "rdp19classified_ORIGINAL.csv.zip"), row.names = 1)
        }

        if (!raw) {
            aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
            counts_original_16 = aligned$counts_original
            counts_original_19 = aligned$counts_original
            matched_taxa_16 <- tax_original16$Taxa[match(colnames(counts_original_16), rownames(tax_original16))]
            matched_taxa_19 <- tax_original19$Taxa[match(colnames(counts_original_19), rownames(tax_original19))]
            colnames(counts_original_16) <- matched_taxa_16
            colnames(counts_original_19) <- matched_taxa_19
            counts_original_16 <- collapse_duplicate_columns_exact(counts_original_16)
            counts_original_19 <- collapse_duplicate_columns_exact(counts_original_19)
            original_names_16 <- colnames(counts_original_16)
            original_names_19 <- colnames(counts_original_19)
            counts_original_16 <- as.data.frame(lapply(counts_original_16, as.numeric), row.names = rownames(counts_original_16), col.names = original_names_16, check.names = FALSE)
            counts_original_19 <- as.data.frame(lapply(counts_original_19, as.numeric), row.names = rownames(counts_original_19), col.names = original_names_19, check.names = FALSE)
            proportions_original_16 <- sweep(counts_original_16, 1, rowSums(counts_original_16), '/')
            proportions_original_19 <- sweep(counts_original_19, 1, rowSums(counts_original_19), '/')
        }
        proportions_original <- sweep(counts_original, 1, rowSums(counts_original), '/')
    } else {
        proportions_original <- NA
        tax_original <- NA
        counts_original_16 <- NA
        counts_original_19 <- NA
        proportions_original_16 <- NA
        proportions_original_19 <- NA
        tax_original16 <- NA
        tax_original19 <- NA
    }

    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {

        temp_dir <- tempfile("repro")
        dir.create(temp_dir)
        unzipped = unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
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
              write.csv(tax_reprocessed2, file = file.path(local, "rdp16classified.csv"), row.names = TRUE)
            } else {
              stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
          }
          
          } else {
            tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
        }

        # ----- Taxonomy reprocessed -----
        unzipped = unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))
        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed= counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
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

        # proportions reprocessed
        proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
        cleanup_tempfiles(temp_dir)
    }

    if (!raw) {
        counts_original = fill_na_zero_numeric(counts_original)
        proportions_original = fill_na_zero_numeric(proportions_original)
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
        proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    }

    if (raw) {
        counts_original_19 <- counts_original
        tax_original19 <- tax_original
        proportions_original_19 <- proportions_original
    }

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = list(rdp19 = counts_original_19, rdp16 = counts_original_16),
            reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
        ),
        proportions = list(
            original = list(rdp19 = proportions_original_19, rdp16 = proportions_original_16),
            reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
        ),
        tax = list(
            original = list(rdp19 = tax_original19, rdp16 = tax_original16),
            reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
        ),
        scale = scale,
        metadata = metadata
    ))
}