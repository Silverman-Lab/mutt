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
    
    # Split scale data into 16S and ITS rows
    scale_16s <- scale %>% 
        select(c("Segment", "Date", "Number", "Sample", "qPCR_log10_16S")) %>%
        mutate(Amplicontype = "16S",
               qPCR_log10_ITS = NA) 
    
    scale_its <- scale %>% 
        select(c("Segment", "Date", "Number", "Sample", "qPCR_log10_ITS")) %>%
        mutate(Amplicontype = "ITS",
               qPCR_log10_16S = NA) 
    
    scale <- bind_rows(scale_16s, scale_its)
    
    metadata = left_join(sra, scale %>% select(c("Segment", "Date", "Number", "Sample", "Amplicontype")), 
                        by = c("Sample", "Amplicontype")) %>%
                        mutate(Sample = paste0(Sample, "_", Amplicontype))
    
    scale = scale %>% 
        select(-c("Segment", "Date", "Number")) %>% 
        mutate(log2_qPCR_16S = qPCR_log10_16S * log2(10)) %>%
        mutate(log2_qPCR_ITS = qPCR_log10_ITS * log2(10)) %>%
        rename(log10_qPCR_16S = qPCR_log10_16S, log10_qPCR_ITS = qPCR_log10_ITS) %>%
        mutate(Sample = paste0(Sample, "_", Amplicontype)) %>% 
        select(-Amplicontype)

    # ----- counts, tax, proportions -----
    counts_16s = read_zipped_table(counts_16s_zip, row.names=NULL) %>%
                    select(-c("Group1","Group2")) %>% rename(Sample = index) %>% 
                    as.data.frame() %>%
                    mutate(Sample = paste0(Sample, "_16S"))
    rownames(counts_16s) = counts_16s$Sample
    counts_16s$Sample = NULL

    counts_ITS = read_zipped_table(counts_ITS_zip, row.names=NULL) %>%
                    select(-c("Group1","Group2")) %>% rename(Sample = index) %>% 
                    as.data.frame() %>%
                    mutate(Sample = paste0(Sample, "_ITS"))
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
    if (file.exists(repro_counts_rds_zip)) {
    temp_dir <- tempfile("repro")
    dir.create(temp_dir)
    unzipped <- unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
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
    unzipped <- unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))
    tax_reprocessed = make_taxa_label(tax_reprocessed)
    
    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
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
                `16S` = counts_16s,
                ITS = counts_ITS
            ),
            reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
        ),
        proportions = list(
            original = list(
                `16S` = proportions_16s,
                ITS = proportions_ITS
            ),
            reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
        ),
        tax = list(
            original = list(
                `16S` = tax_16s,
                ITS = tax_ITS
            ),
            reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
        ),
        scale = scale,
        metadata = metadata
    ))


}