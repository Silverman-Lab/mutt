parse_2024_jin_pnas_semen <- function(raw = FALSE, align = FALSE) {
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

  library(tidyverse)
  library(readxl)
  library(stringr)
  library(readr)
  
  # -----local path ---------------------------
  local <- file.path("2024_jin_pnas_semen")

  # ---------- file paths ---------------------
  repro_counts_rds_zip    <- file.path(local, "PRJNA747100_dada2_counts.rds.zip")
  repro_tax_zip           <- file.path(local, "PRJNA747100_dada2_taxa.rds.zip")
  metadata_zip            <- file.path(local, "metadata4.txt.zip")
  sra_zip                 <- file.path(local, "SraRunTable (38).csv.zip")
  cfu_zip                 <- file.path(local, "CFU_data.txt.zip")
  counts_zip              <- file.path(local, "table.tsv.zip")
  tax_zip                 <- file.path(local, "taxonomy.tsv.zip")
  
  # ---- metadata ----
  metadata = read_zipped_table(metadata_zip, sep="\t", row.names=NULL) %>% as.data.frame() %>% rename(Sample = `sample-ID`)
  metadata <- metadata[!(metadata$Sample %in% c("S1", "S13", "SN", "SNTC")), ]
  sra = read_zipped_table(sra_zip, row.names=NULL) %>% as.data.frame() 
  sra = sra %>% rename(Sample = `Sample Name`, Accession = Run)
  metadata = full_join(metadata, sra, by = "Sample")

  # ---- scale ----
  scale <- read_zipped_table(cfu_zip, sep="\t", row.names=NULL) %>% as.data.frame() %>% rename(Sample = `sample-ID`)
  scale <- scale %>% mutate(log2_CFU = ifelse(CFU > 0, log2(CFU),NA)) %>% 
                  mutate(log10_CFU = ifelse(CFU > 0, log10(CFU),NA))
  
  # ---- counts ----
  counts = read_zipped_table(counts_zip, sep="\t", row.names=NULL) %>% as.data.frame() %>% rename(OTUID = `OTU ID`)
  rownames(counts) <- counts$OTUID
  counts <- subset(counts, select = -OTUID)
  counts <- t(counts)
  counts <- counts[!(rownames(counts) %in% c("S1", "SN", "SNTC")), ] %>% as.data.frame()
  
  # ---- tax ----
  tax <- read_zipped_table(tax_zip, sep="\t", row.names=NULL) %>% as.data.frame() %>%
          rename(OTUID = `Feature ID`)

  tax <- tax %>%
  separate(Taxon,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";", fill = "right") %>%
  mutate(across(Kingdom:Species, ~ trimws(.))) %>% 
  mutate(across(Kingdom:Species, ~ gsub("^(k__|p__|c__|o__|f__|g__|s__)", "", .))) %>%
  mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), "unclassified", .)))
  tax <- make_taxa_label(tax)
  rownames(tax) <- tax$OTUID
  if (!raw) {
      matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
      colnames(counts) <- matched_taxa
      counts <- collapse_duplicate_columns_exact(counts)
      original_names <- colnames(counts)
      counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
  }
  # --- proportions ---
  proportions_original <- sweep(counts, MARGIN = 1,STATS  = rowSums(counts), FUN = "/")
  
  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
    temp_dir <- tempfile("repro")
    dir.create(temp_dir)
    unzipped = unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
    counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(counts_file))

    # ----- Taxonomy reprocessed -----
    unzipped = unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))
    
    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
    aligned = rename_and_align(counts_reprocessed = counts_reprocessed, counts_original = counts, 
                              proportions_original = proportions_original, metadata=metadata, scale=scale, 
                              by_col="Sample", align = align, study_name=basename(local))
      counts_reprocessed = aligned$reprocessed
      counts_original = aligned$counts_original
      proportions_original = aligned$proportions_original
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
      original_names <- colnames(counts_reprocessed)
      counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
      original_names <- colnames(counts_original)
      counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
      original_names <- colnames(proportions_original)
      proportions_original <- as.data.frame(lapply(proportions_original, as.numeric), row.names = rownames(proportions_original), col.names = original_names, check.names = FALSE)
    }

    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
    cleanup_tempfiles(temp_dir)
  }

  if (!raw) {
      counts = fill_na_zero_numeric(counts)
      proportions_original = fill_na_zero_numeric(proportions_original)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }
  
  return(list(scale=scale, 
              metadata=metadata, 
              counts=list(
                original = counts, 
                reprocessed = counts_reprocessed
              ),
              proportions=list(
                original = proportions_original,
                reprocessed = proportions_reprocessed
              ),
              tax=list(
                original = tax,
                reprocessed = tax_reprocessed
              )
  )
  )
}