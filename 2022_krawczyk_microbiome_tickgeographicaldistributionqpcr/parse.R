parse_2022_krawczyk_microbiome_tickgeographicaldistributionqpcr <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  if (!is.logical(raw) || length(raw) != 1) {
    stop("raw must be a boolean")
  }
  if (!is.logical(align) || length(align) != 1) {
    stop("align must be a boolean")
  }
  
  library(tibble)
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2022_krawczyk_microbiome_tickgeographicaldistributionqpcr")

  # ----- File paths -----
  metadata_two_zip     <- file.path(local, "40168_2022_1276_MOESM1_ESM.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  counts_zip           <- file.path(local, "originalcounts.csv.zip")
  scale_zip            <- file.path(local, "scale_qpcr.csv.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA813158_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA813158_dada2_taxa.rds.zip")

  # ----- Scale -----
  scale = read_zipped_table(scale_zip, row.names = NULL) %>% as.data.frame() %>% mutate(qpcr_16s_ng_ul = as.numeric(`16S rRNA content in ng/µL`)) %>%
    mutate(log2_qPCR_16S_ng_ul = ifelse(qpcr_16s_ng_ul > 0, log2(qpcr_16s_ng_ul), NA)) %>%
    mutate(log10_qPCR_16S_ng_ul = ifelse(qpcr_16s_ng_ul > 0, log10(qpcr_16s_ng_ul), NA)) %>% 
    rename(Sample_name = `Sample ID`) %>% select(Sample_name, qpcr_16s_ng_ul, log2_qPCR_16S_ng_ul, log10_qPCR_16S_ng_ul)

  # ----- Metadata -----
  metadata  = read_zipped_table(metadata_zip, row.names = NULL) %>% as.data.frame() %>% rename(Sample_name = `Sample ID`)
  meta_two_file <- unzip(metadata_two_zip, list = TRUE)$Name[1]
  extracted_xlsx <- unzip(metadata_two_zip, files = meta_two_file, exdir = tempdir(), overwrite = TRUE)[1]
  metadata_two <- readxl::read_xlsx(extracted_xlsx, sheet = 2, col_names = TRUE) %>% rename(Sample_name = `Sample ID`)
  metadata <- left_join(metadata, metadata_two, by = "Sample_name") %>% 
              rename(Accession = Run)

  # -------- Counts --------
  if (file.exists(counts_zip)) {
    counts_original = read_zipped_table(counts_zip, row.names = NULL)
    columns_to_drop <- c("Name","Taxonomy", "Combined Abundance",	"Min",	"Max",	"Mean",	"Median",	"Std")

    # Taxa
    if (!file.exists(file.path(local,"taxawsilvaandrdp.csv.zip"))) {
      taxonomy_cols <- c("Sequence","Taxonomy")
      tax <- counts_original[, names(counts_original) %in% taxonomy_cols]
      tax <- tax %>%
        separate(
          col = Taxonomy,
          into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
          sep = ",\\s*",
          remove = FALSE, 
          extra = "merge",   # merge extra pieces into Species
          fill = "right"    # fill missing with NA
        ) %>%
        mutate(across(Kingdom:Species, ~ sub("D_\\d+__", "", .)))
      tax$species <- tax$Species
      tax$Species <- NULL
      tax = make_taxa_label(tax) 
      tax$Species <- tax$species
      tax$species <- NULL
      if (file.exists(file.path("helperdata/rdp_19_toGenus_trainset.fa.gz"))) {
        required_pkgs <- c("dada2", "Biostrings")
        missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
        if (length(missing_pkgs) > 0) {
          stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
              ". Please install them before running this function.")
        }
        tax = tax %>% rename_with(~ paste0(., "_silva"), .cols = c(Taxa, Kingdom, Phylum, Class, Order, Family, Genus, Species))
        seqs <- Biostrings::DNAStringSet(tax$Sequence)
        rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_19_toGenus_trainset.fa.gz"), multithread=TRUE) %>% as.data.frame()
        tax <- cbind(tax, rdpclassified)
        tax = make_taxa_label(tax) 
      }
      # Counts
      rownames(tax) <- tax$Sequence
      write.csv(tax, file = file.path(local, "taxawsilvaandrdp.csv"), row.names = TRUE)
    } else {
      tax <- read_zipped_table(file.path(local, "taxawsilvaandrdp.csv.zip"), row.names = 1)
    }

    counts_original <- counts_original %>% column_to_rownames("Sequence")
    counts_original = counts_original[, !(names(counts_original) %in% columns_to_drop)]
    counts_original = as.data.frame(t(counts_original))
    
    if (!raw) {
      aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
      counts_original = aligned$counts_original
      matched_taxa <- tax$Taxa[match(colnames(counts_original), rownames(tax))]
      colnames(counts_original) <- matched_taxa
      counts_original <- collapse_duplicate_columns_exact(counts_original)
      original_names <- colnames(counts_original)
      counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
    }
    # Calculate proportions
    proportions_original <- sweep(counts_original, MARGIN = 1, STATS = rowSums(counts_original), FUN = "/")
  } else {
    counts_original = NA
    proportions_original = NA
  }

    if (file.exists(counts_zip)) {
    counts_original2 = read_zipped_table(counts_zip, row.names = NULL)
    columns_to_drop <- c("Name","Taxonomy", "Combined Abundance",	"Min",	"Max",	"Mean",	"Median",	"Std")

    # Taxa
    if (!file.exists(file.path(local,"taxawsilvaandrdp2.csv.zip"))) {
      taxonomy_cols <- c("Sequence","Taxonomy")
      tax2 <- counts_original2[, names(counts_original2) %in% taxonomy_cols]
      tax2 <- tax2 %>%
        separate(
          col = Taxonomy,
          into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
          sep = ",\\s*",
          remove = FALSE, 
          extra = "merge",   # merge extra pieces into Species
          fill = "right"    # fill missing with NA
        ) %>%
        mutate(across(Kingdom:Species, ~ sub("D_\\d+__", "", .)))
      tax2$species <- tax2$Species
      tax2$Species <- NULL
      tax2 = make_taxa_label(tax2) 
      tax2$Species <- tax2$species
      tax2$species <- NULL
      if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
        required_pkgs <- c("dada2", "Biostrings")
        missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
        if (length(missing_pkgs) > 0) {
          stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
              ". Please install them before running this function.")
        }
        tax2 = tax2 %>% rename_with(~ paste0(., "_silva"), .cols = c(Taxa, Kingdom, Phylum, Class, Order, Family, Genus, Species))
        seqs <- Biostrings::DNAStringSet(tax2$Sequence)
        rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
        tax2 <- cbind(tax2, rdpclassified)
        tax2 = make_taxa_label(tax2) 
      }
      # Counts
      rownames(tax2) <- tax2$Sequence
      write.csv(tax2, file = file.path(local, "taxawsilvaandrdp2.csv"), row.names = TRUE)
    } else {
      tax2 <- read_zipped_table(file.path(local, "taxawsilvaandrdp2.csv.zip"), row.names = 1)
    }

    counts_original2 <- counts_original2 %>% column_to_rownames("Sequence")
    counts_original2 = counts_original2[, !(names(counts_original2) %in% columns_to_drop)]
    counts_original2 = as.data.frame(t(counts_original2))
    
    if (!raw) {
      aligned = rename_and_align(counts_original = counts_original2, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
      counts_original2 = aligned$counts_original2
      matched_taxa <- tax2$Taxa[match(colnames(counts_original2), rownames(tax2))]
      colnames(counts_original2) <- matched_taxa
      counts_original2 <- collapse_duplicate_columns_exact(counts_original2)
      original_names <- colnames(counts_original2)
      counts_original2 <- as.data.frame(lapply(counts_original2, as.numeric), row.names = rownames(counts_original2), col.names = original_names, check.names = FALSE)
    }
    # Calculate proportions
    proportions_original2 <- sweep(counts_original2, MARGIN = 1, STATS = rowSums(counts_original2), FUN = "/")
  } else {
    counts_original2 = NA
    proportions_original2 = NA
  }

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
        NULL
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
    }

    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
    proportions_reprocessed2 <- sweep(counts_reprocessed2, 1, rowSums(counts_reprocessed2), '/')


    cleanup_tempfiles(temp_dir)

    if (!raw) {
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
        proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    }
  } else {
    counts_reprocessed = NA
    proportions_reprocessed = NA
    tax_reprocessed = NA
    counts_reprocessed2 = NA
    proportions_reprocessed2 = NA
    tax_reprocessed2 = NA
  }



  # ----- Return all -----
  return(list(
    counts      = list(
      original = list(
                rdp19 = counts_original,
                rdp16 = counts_original2
      ),
      reprocessed = list(
                rdp19 = counts_reprocessed,
                rdp16 = counts_reprocessed2
      )
    ),
    proportions = list(
      original = list(
                rdp19 = proportions_original,
                rdp16 = proportions_original2
      ),
      reprocessed = list(
                rdp19 = proportions_reprocessed,
                rdp16 = proportions_reprocessed2
      )
    ),
    tax         = list(
      original = list(
                rdp19 = tax,
                rdp16 = tax2
      ),
      reprocessed = list(
                rdp19 = tax_reprocessed,
                rdp16 = tax_reprocessed2
      )
    ),
    scale       = scale,
    metadata    = metadata
  ))
}
