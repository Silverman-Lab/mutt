<<<<<<< Updated upstream
parse_2022_krawczyk_microbiome_tickgeographicaldistributionqpcr <- function(raw = FALSE, align = FALSE) {
=======
parse_2022_krawczyk_microbiome_tickgeographicaldistributionqpcr <- function() {
>>>>>>> Stashed changes
  required_pkgs <- c("tibble", "tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
<<<<<<< Updated upstream
  if (!is.logical(raw) || length(raw) != 1) {
    stop("raw must be a boolean")
  }
  if (!is.logical(align) || length(align) != 1) {
    stop("align must be a boolean")
  }
  
=======
>>>>>>> Stashed changes
  library(tibble)
  library(tidyverse)
  library(readxl)

<<<<<<< Updated upstream
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
      counts_original2 = aligned$counts_original
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
=======
  # ----- File paths -----
  local               <- "2022_krawczyk_microbiome_tickgeographicaldistributionqpcr/"
  metadata_two_zip    <- paste0(local, "40168_2022_1276_MOESM1_ESM.zip")
  metadata_zip        <- paste0(local, "SraRunTable (30).csv")
  counts_zip          <- paste0(local, "40168_2022_1276_MOESM3_ESM.zip")
  scale_zip           <- paste0(local, "scale_qpcr.zip")
  repro_counts_rds_zip<-  paste0(local, "PRJNA813158_dada2_merged_nochim.rds.zip"
  repro_tax_zip       <-  paste0(local, "PRJNA813158_dada2_taxonomy_merged.rds.zip"

  # ----- Scale -----
  if (file.exists(scale_zip)) {
    scale_file <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_con  <- unz(scale_zip, scale_file)
    scaledata  <- read.csv(scale_con) %>% as.data.frame()
    cleaned <- as.numeric(scaledata$`16S rRNA content in ng/µL`)  # Will coerce character "no DNA left" to NA
    scale = data.frame(Sample = scaledata$`Sample ID`, `16S rRNA content in ng/µL` = cleaned)
  } else {
    scale = NA
  }
  # ----- Metadata -----
  if (file.exists(metadata_zip)) {
    meta_file <- unzip(metadata_zip, list = TRUE)$Name[1]
    meta_con  <- unz(metadata_zip, meta_file)
    metadata  = read.csv(meta_con) %>% as.data.frame()
    meta_two_file <- unzip(metadata_two_zip, list = TRUE)$Name[1]
    extracted_xlsx <- unzip(metadata_two_zip, files = meta_two_file, exdir = tempdir(), overwrite = TRUE)[1]
    metadata_two <- read_excel(extracted_xlsx, sheet = 2)
    metadata <- left_join(metadata, metadata_two, by = "Sample ID")
  } else {
    metadata = NA
  }

  # -------- Counts --------
  if (file.exists(counts_zip)) {
    counts_file <- unzip(counts_zip, list = TRUE)$Name[1]
    extracted_xlsx <- unzip(counts_zip, files = counts_file, exdir = tempdir(), overwrite = TRUE)[1]
    counts_original = read_excel(extracted_xlsx, sheet = 1)
    columns_to_drop <- c("Name", "Taxonomy", "Combined Abundance",	"Min,"	"Max",	"Mean",	"Median",	"Std")
    
    # Taxa
    taxonomy_cols <- c("Name", "Taxonomy")
    tax <- counts_original[, names(counts_original) %in% taxonomy_cols]
    rownames(tax) <- paste0("Taxon_", seq_len(nrow(tax)))
    tax <- tax %>%
      separate(
        col = Taxonomy,
        into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
        sep = ",\\s*",
        remove = FALSE
      ) %>%
      mutate(across(Kingdom:Species, ~ sub("D_\\d+__", "", .)))

    # Counts
    counts_original = counts_original[, !(names(counts_original) %in% columns_to_drop)]
    rownames(counts_original) <- paste0("Taxon_", seq_len(nrow(counts_original)))
  } else {
    counts_original = NA
  }

  if (file.exists(orig_counts_zip)) {
    proportions_original = counts_original
    proportions_original[-1] <- lapply(
      counts_original[-1],
      function(col) col / sum(col)
    )
  } else if (file.exists(orig_prop_zip)) {
    prop_csv <- unzip(orig_prop_zip, list = TRUE)$Name[1]
    prop_con <- unz(orig_prop_zip, prop_csv)
    proportions_original = read.csv(prop_con, row.names = 1, check.names = FALSE) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sequence") %>%
      dplyr::select(Sequence, everything())
  } else {
  proportions_original <- NA
  }

  # ----- Reprocessed counts from RDS ZIP -----
  temp_rds            <- tempfile(fileext = ".rds")
  unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)
  rds_file            <- list.files(dirname(temp_rds), pattern = "\\.rds$", full.names = TRUE)[1]
  seqtab_nochim       <- readRDS(rds_file)
  rpt_mat             <- t(seqtab_nochim)
  counts_reprocessed  <- as.data.frame(rpt_mat)
  counts_reprocessed$Sequence <- rownames(counts_reprocessed)
  counts_reprocessed = counts_reprocessed[, c("Sequence", setdiff(names(counts_reprocessed), "Sequence"))]
  rownames(counts_reprocessed) <- paste0("Taxon_", seq_len(nrow(counts_reprocessed)))

  # proportions reprocessed
  proportions_reprocessed = counts_reprocessed
  proportions_reprocessed[-1] <- lapply(
    counts_reprocessed[-1],
    function(col) col / sum(col)
  )

  # ----- Taxonomy reprocessed -----
  temp_tax <- tempfile(fileext = ".rds")
  unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
  tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
  taxonomy_matrix <- readRDS(tax_file)
  rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
  tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
  tax_reprocessed = tax_table

  # ----- Return all -----
  return(list(
    counts      = list(
      original = counts_original
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original = proportions_original
      reprocessed = proportions_reprocessed
    ),
    tax         = list(
      original = tax_original
      reprocessed = tax_reprocessed
    )
>>>>>>> Stashed changes
    scale       = scale,
    metadata    = metadata
  ))
}
