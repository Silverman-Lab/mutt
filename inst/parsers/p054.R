parse_2024_sternes_frontmicrobiol_IBDppiqPCR <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl")
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
  
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2024_sternes_frontmicrobiol_IBDppiqPCR")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "PRJNA1120972_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1120972_dada2_taxa.rds.zip")
  metadata_zip         <- file.path(local, "SraRunTable (34).csv.zip")
  scale_zip            <- file.path(local, "scale.csv.zip")

  # --- original counts, proportions, tax ---
  counts <- NA
  proportions <- NA
  tax <- NA
  
  # --- metadata and scale ----
  metadata     <- read_zipped_table(metadata_zip, row.names=NULL) %>% as.data.frame() %>% 
                  rename(Sample_name = sampleid, Accession = Run) %>%
                  mutate(Sample_name = as.character(Sample_name))
  scale_raw    <- read_zipped_table(scale_zip, row.names=NULL) %>% as.data.frame() %>% 
                  rename(Sample_name = ID) %>%  mutate(Sample_name = as.character(Sample_name)) %>% as.data.frame() %>% 
                  mutate(log2_DU_bacterial_load = ifelse(DU_bacterial_load > 0, log2(DU_bacterial_load),NA),
                         log2_TI_bacterial_load = ifelse(TI_bacterial_load > 0, log2(TI_bacterial_load),NA),
                         log2_DNS_bacterial_load = ifelse(DNS_bacterial_load > 0, log2(DNS_bacterial_load),NA)) %>% 
                  mutate(log10_DU_bacterial_load = ifelse(DU_bacterial_load > 0, log10(DU_bacterial_load),NA),
                         log10_TI_bacterial_load = ifelse(TI_bacterial_load > 0, log10(TI_bacterial_load),NA),
                         log10_DNS_bacterial_load = ifelse(DNS_bacterial_load > 0, log10(DNS_bacterial_load),NA))

  meta_part    <- scale_raw %>% select(c("Sample_name", "DX_Groups","AGE","age_binned","Sex", 
  "NC_total_score", "NC_binned", "BMI",	"BMI_binned", "SAGIS_binned", "SAGIS_binned2", "PPI")) 
  scale        = scale_raw %>% select(-c("DX_Groups","AGE","age_binned","Sex", 
  "NC_total_score", "NC_binned", "BMI",	"BMI_binned", "SAGIS_binned", "SAGIS_binned2", "PPI")) 
  metadata     = metadata %>% left_join(meta_part, by = "Sample_name")

  metadata = metadata %>%
  mutate(
    across(
      c(age_binned, NC_binned, BMI_binned, SAGIS_binned, SAGIS_binned2, Sex, PPI),
      ~ na_if(.x, "")
    )
  ) %>%
  
  mutate(
    NC_binned    = recode(NC_binned, `1000-2000` = "1000+"),
    SAGIS_binned = recode(SAGIS_binned, `10/19/25` = "10-19"),
    BMI_binned   = recode(BMI_binned, Severely_Obese = "Severely Obese"),
    DX_Groups    = recode(DX_Groups,
                          KIT18CONT = "CONTROL",
                          KIT19CONT = "CONTROL")
  ) %>%
  
  mutate(
    age_binned   = factor(age_binned,
                          levels = c("17-29","30-39","40-49","50-59","60-69","70-79"),
                          ordered = TRUE),
    NC_binned    = factor(NC_binned,
                          levels = c("0-499","500-999","1000+"),
                          ordered = TRUE),
    BMI_binned   = factor(BMI_binned,
                          levels = c("Underweight","Normal","Overweight","Obese","Severely Obese"),
                          ordered = TRUE),
    SAGIS_binned = factor(SAGIS_binned,
                          levels = c("0-9","10-19","20-29","30-80"),
                          ordered = TRUE),
    SAGIS_binned2= factor(SAGIS_binned2,
                          levels = c("Low","High"),
                          ordered = TRUE),
    Sex          = factor(Sex, levels = c("M","F")),
    PPI          = factor(PPI, levels = c("No PPI","PPI")),
    DX_Groups    = factor(DX_Groups)
  )

  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA 
  tax_reprocessed2 <- NA
  counts_reprocessed2 <- NA
  proportions_reprocessed2 <- NA

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
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

    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name",align = align, study_name=basename(local))
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
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
      proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
  }
  
  return(list(
    counts = list(
      original = NA,
      reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
    ),
    proportions = list(
      original = NA,
      reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
    ),
    tax = list(
      original = NA,
      reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
    ),
    scale = scale,
    metadata = metadata
  ))
}