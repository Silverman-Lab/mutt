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
  dataset <- read_zipped_table(counts_zip, row.names = NULL)

  # ---- metadata ----
  metadata_cols <- c("Subject", "Timepoint", "Milling", "Frequency", "StoolsperDay",
                     "BristolStoolScale_highest", "WaterContent_perc", "pH",
                     "Calprotectin_ugperg", "MPO_ngperml")
  metadata <- dataset[, metadata_cols]
  sra = read_zipped_table(metadata_zip, row.names=NULL) %>% rename(Accession = Run)
  metadata1 <- read_zipped_table(metadata_zip_1, row.names = NULL) %>%
    mutate(Replicate = gsub("Replicate ", "Replicate_", ID))
  metadata2 <- read_zipped_table(metadata_zip_2, row.names = NULL) %>%
    pivot_longer(cols = starts_with("Replicate"), names_to = "Replicate", values_to = "Value") %>%
    pivot_wider(names_from = ID, values_from = Value)

  merged_metadata <- metadata1 %>%
    full_join(metadata2, by = "Replicate") #%>% full_join(metadata, by = c("Subject", "Timepoint"))

  # ---- scale ----
  scale_cols <- c("Subject", "Timepoint",
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
  counts <- bind_cols(dataset[, c("Subject", "Timepoint")], counts_original)
  counts_original <- counts %>% select(-c(3:19))
  counts_metabolomics <- counts %>% select(c(1:19))

  # ---- tax ----
  tax <- tibble(
    taxonomy = colnames(counts_original)[!(colnames(counts_original) %in% c(metadata_cols, scale_cols))]
  )

  # ---- proportions ----
  row_sums <- rowSums(counts_original)
  prop_mat <- sweep(as.matrix(counts_original), 1, row_sums, FUN = "/")
  prop_mat[is.nan(prop_mat)] <- 0
  proportions <- bind_cols(dataset[, c("Subject", "Timepoint")], as_tibble(prop_mat))

  row_sums <- rowSums(counts_metabolomics)
  prop_mat <- sweep(as.matrix(counts_metabolomics), 1, row_sums, FUN = "/")
  prop_mat[is.nan(prop_mat)] <- 0
  proportions_metabolomics <- bind_cols(dataset[, c("Subject", "Timepoint")], as_tibble(prop_mat))

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
  aligned = rename_and_align(counts_reprocessed = counts_reprocessed, counts_original = counts_original, 
                            proportions_original = proportions, metadata=metadata, scale=scale, 
                            by_col=c("Subject", "Timepoint"), align = align, study_name=basename(local))
  counts_reprocessed = aligned$reprocessed
  counts_original = aligned$counts_original
  proportions = aligned$proportions_original
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
    tax         = list(original = list(metabolomics = NA, amplicon = tax), reprocessed = list(metabolomics = NA, amplicon = tax_reprocessed)),
    scale       = scale,
    metadata    = list(metadata1 = merged_metadata, metadata2 = metadata, sra=sra)
  ))
}