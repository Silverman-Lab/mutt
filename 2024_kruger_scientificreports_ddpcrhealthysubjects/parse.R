parse_2024_kruger_scientificreports_ddpcrhealthysubjects <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tibble", "tidyverse")
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

  library(tibble)
  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2024_kruger_scientificreports_ddpcrhealthysubjects")

  # ----- File paths -----
  counts_zip           <- file.path(local, "41598_2024_75477_MOESM2_ESM.csv.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  metadata_zip_1       <- file.path(local, "41598_2024_75477_MOESM3_ESM.csv.zip")
  metadata_zip_2       <- file.path(local, "41598_2024_75477_MOESM4_ESM.csv.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA1162476_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1162476_dada2_taxa.rds.zip")


  # ---- counts ----
  if (!file.exists(counts_zip)) stop("Counts file not found: ", counts_zip)
  dataset <- read_zipped_table(counts_zip, row.names = NULL) %>% 
    mutate(Frequency = gsub("Day ", "", Frequency)) %>%
    mutate(subject_timepoint_replicate = ifelse(is.na(Frequency) | Frequency == "", 
                                              paste(Subject, Timepoint, sep = "_"),
                                              paste(Subject, Timepoint, Frequency, sep = "_")))
  colnames(dataset) <- ifelse(
    grepl("^[kpcofg]__[^_]*$", colnames(dataset)),  
    gsub("([kpcofg])__", "\\1_", colnames(dataset)),  
    colnames(dataset)  
  )

  # ---- metadata ----
  metadata_cols <- c("subject_timepoint_replicate", "Subject", "Timepoint", "Milling", "Frequency", "StoolsperDay",
                     "BristolStoolScale_highest", "WaterContent_perc", "pH",
                     "Calprotectin_ugperg", "MPO_ngperml", "PhylogeneticDiversity", "Chao1", "inverse_simpson", "gini_simpson", "shannon", "fisher")
  metadata <- dataset[, metadata_cols] 
  sra = read_zipped_table(metadata_zip, row.names=NULL) %>% 
    rename(Accession = Run, Sample = `Sample Name`) %>%
    mutate(
      Subject = paste0(gsub(".*Volunteer(\\d+).*", "\\1", `description:_replicate`)),
      Timepoint = gsub("^S\\d+\\.(\\d+).*", "\\1", Sample),
      Replicate = ifelse(grepl("Replicate", `description:_replicate`), Timepoint, NA),
      Processing = case_when(
        grepl("mill homogenized", `description:_sample_information`) ~ "mill_homogenized",
        grepl("hammered not milled", `description:_sample_information`) ~ "hammered",
        TRUE ~ NA_character_
      ),
      subject_timepoint_replicate = ifelse(is.na(Replicate) | Replicate == "", 
                                              paste(Subject, Timepoint, sep = "_"),
                                              paste(Subject, Timepoint, Replicate, sep = "_"))
    )
  metadata1 <- read_zipped_table(metadata_zip_1, row.names = NULL) %>%
    mutate(Replicate = gsub("Replicate ", "Replicate_", ID))
  metadata2 <- read_zipped_table(metadata_zip_2, row.names = NULL) %>%
    pivot_longer(cols = starts_with("Replicate"), names_to = "Replicate", values_to = "Value") %>%
    pivot_wider(names_from = ID, values_from = Value)

  merged_metadata <- metadata1 %>%
    full_join(metadata2, by = "Replicate") 

  # ---- scale ----
  scale_cols <- c("subject_timepoint_replicate",
                  "Mean Fungi copies_per mg total weight", 
                  "Mean Fungi copies_per mg dry weight", 
                  "Mean bacteria copies_per mg total weight", 
                  "Mean bacteria copies_per mg dry weight")
  scale    <- dataset[, scale_cols] 
  scale <- scale %>% mutate(log2_Mean_Fungi_copies_per_mg_total_weight = ifelse(`Mean Fungi copies_per mg total weight` > 0, log2(`Mean Fungi copies_per mg total weight`),NA),
                         log2_Mean_Fungi_copies_per_mg_dry_weight = ifelse(`Mean Fungi copies_per mg dry weight` > 0, log2(`Mean Fungi copies_per mg dry weight`),NA),
                         log2_Mean_bacteria_copies_per_mg_total_weight = ifelse(`Mean bacteria copies_per mg total weight` > 0, log2(`Mean bacteria copies_per mg total weight`),NA),
                         log2_Mean_bacteria_copies_per_mg_dry_weight = ifelse(`Mean bacteria copies_per mg dry weight` > 0, log2(`Mean bacteria copies_per mg dry weight`),NA),
                         log10_Mean_Fungi_copies_per_mg_total_weight = ifelse(`Mean Fungi copies_per mg total weight` > 0, log10(`Mean Fungi copies_per mg total weight`),NA),
                         log10_Mean_Fungi_copies_per_mg_dry_weight = ifelse(`Mean Fungi copies_per mg dry weight` > 0, log10(`Mean Fungi copies_per mg dry weight`),NA),
                         log10_Mean_bacteria_copies_per_mg_total_weight = ifelse(`Mean bacteria copies_per mg total weight` > 0, log10(`Mean bacteria copies_per mg total weight`),NA),
                         log10_Mean_bacteria_copies_per_mg_dry_weight = ifelse(`Mean bacteria copies_per mg dry weight` > 0, log10(`Mean bacteria copies_per mg dry weight`),NA)) 

  scale = scale %>% left_join(sra %>% select(Sample, subject_timepoint_replicate), by = "subject_timepoint_replicate")
  
  scale_cols <- c("subject_timepoint_replicate", "TotalSCFAs_umolpermg", "TotalBCFAs_umolpermg" )
  metabolomicsscale <- dataset %>% 
    select(scale_cols) %>%
    mutate(log2_TotalSCFAs_umolpermg = ifelse(`TotalSCFAs_umolpermg` > 0, log2(`TotalSCFAs_umolpermg`),NA),
           log2_TotalBCFAs_umolpermg = ifelse(`TotalBCFAs_umolpermg` > 0, log2(`TotalBCFAs_umolpermg`),NA),
           log10_TotalSCFAs_umolpermg = ifelse(`TotalSCFAs_umolpermg` > 0, log10(`TotalSCFAs_umolpermg`),NA),
           log10_TotalBCFAs_umolpermg = ifelse(`TotalBCFAs_umolpermg` > 0, log10(`TotalBCFAs_umolpermg`),NA))

  counts_original <- dataset %>% 
    select(-c(3:31)) %>%
    column_to_rownames("subject_timepoint_replicate") %>%
    select(-c("Subject", "Timepoint"))
    
  counts_metabolomics <- dataset %>% 
    select(subject_timepoint_replicate, c(13:23)) %>%
    column_to_rownames("subject_timepoint_replicate")

  # ---- tax ----
  tax <- tibble(
    taxonomy = colnames(counts_original)[!(colnames(counts_original) %in% c(metadata_cols, scale_cols))]
  )

  tax_metabolomics <- tibble(
    metabolites = colnames(counts_metabolomics)[!(colnames(counts_metabolomics) %in% c(metadata_cols, scale_cols))]
  )

  # ---- proportions ----
  proportions <- sweep(counts_original, 1, rowSums(counts_original), '/')
  proportions_metabolomics <- sweep(counts_metabolomics, 1, rowSums(counts_metabolomics), '/')

  if (!raw) {
    aligned = rename_and_align(counts_original = counts_original, 
                              proportions_original = proportions, metadata=metadata, scale=scale, 
                              by_col="subject_timepoint_replicate", align = align, study_name=basename(local))
    counts_original = aligned$counts_original
    proportions = aligned$proportions_original
    original_names <- colnames(counts_original)
    counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
    original_names <- colnames(proportions)
    proportions <- as.data.frame(lapply(proportions, as.numeric), row.names = rownames(proportions), col.names = original_names, check.names = FALSE)
    aligned = rename_and_align(counts_original = counts_metabolomics, 
                              proportions_original = proportions_metabolomics, metadata=metadata, scale=scale, 
                              by_col="subject_timepoint_replicate", align = align, study_name=basename(local))
    counts_metabolomics = aligned$counts_original
    proportions_metabolomics = aligned$proportions_original
    original_names <- colnames(counts_metabolomics)
    counts_metabolomics <- as.data.frame(lapply(counts_metabolomics, as.numeric), row.names = rownames(counts_metabolomics), col.names = original_names, check.names = FALSE)
    original_names <- colnames(proportions_metabolomics)
    proportions_metabolomics <- as.data.frame(lapply(proportions_metabolomics, as.numeric), row.names = rownames(proportions_metabolomics), col.names = original_names, check.names = FALSE)
  }

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
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
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=sra, scale=scale, by_col="Sample", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
        original_names <- colnames(counts_reprocessed)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
    }

    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
    cleanup_tempfiles(temp_dir)
  }

  if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      counts_metabolomics = fill_na_zero_numeric(counts_metabolomics)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions = fill_na_zero_numeric(proportions)
      proportions_metabolomics = fill_na_zero_numeric(proportions_metabolomics)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

  return(list(
    counts      = list(original = list(metabolomics = counts_metabolomics, amplicon = counts_original), reprocessed = list(metabolomics = NA, amplicon = counts_reprocessed)),
    proportions = list(original = list(metabolomics = proportions_metabolomics, amplicon = proportions), reprocessed = list(metabolomics = NA, amplicon = proportions_reprocessed)),
    tax         = list(original = list(metabolomics = tax_metabolomics, amplicon = tax), reprocessed = list(metabolomics = NA, amplicon = tax_reprocessed)),
    scale       = scale,
    metadata    = list(replicates = merged_metadata, dataset = metadata, sra = sra)
  ))
}