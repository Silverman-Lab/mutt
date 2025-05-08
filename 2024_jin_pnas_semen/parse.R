parse_2024_jin_pnas_semen <- function(raw=FALSE, align=FALSE) {
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
      counts <- as.data.frame(t(rowsum(t(counts), group = colnames(counts))))
  }
  # --- proportions ---
  proportions_original <- sweep(counts, MARGIN = 1,STATS  = rowSums(counts), FUN = "/")
  
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
  }

  # taxa
  if (!raw) {
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
  }

  # proportions reprocessed
  proportions_reprocessed = counts_reprocessed
  proportions_reprocessed[-1] <- lapply(
      counts_reprocessed[-1],
      function(col) col / sum(col)
  )

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