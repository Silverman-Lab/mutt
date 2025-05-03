parse_2024_kruger_scientificreports_ddpcrhealthysubjects <- function() {
  required_pkgs <- c("tibble", "tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)

<<<<<<< Updated upstream
<<<<<<< Updated upstream
  # ----- Local base directory -----
  local <- file.path("2024_kruger_scientificreports_ddpcrhealthysubjects")

  # ----- File paths -----
  counts_zip           <- file.path(local, "41598_2024_75477_MOESM2_ESM.csv.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  metadata_zip_1       <- file.path(local, "41598_2024_75477_MOESM3_ESM.csv.zip")
  metadata_zip_2       <- file.path(local, "41598_2024_75477_MOESM4_ESM.csv.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA1162476_dada2_merged_nochim.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1162476_dada2_taxonomy_merged.rds.zip")

=======
  local               <- "2024_kruger_scientificreports_ddpcrhealthysubjects/"
  counts_zip          <- paste0(local, "41598_2024_75477_MOESM2_ESM.csv.zip")
  metadata_zip        <- paste0(local, "SraRunTable (30).csv.zip")
  metadata_zip_1      <- paste0(local, "41598_2024_75477_MOESM3_ESM.csv.zip")
  metadata_zip_2      <- paste0(local, "41598_2024_75477_MOESM4_ESM.csv.zip")
  repro_counts_rds_zip<- paste0(local, "PRJNA1162476_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- paste0(local, "PRJNA1162476_dada2_taxonomy_merged.rds.zip")

>>>>>>> Stashed changes
=======
  local               <- "2024_kruger_scientificreports_ddpcrhealthysubjects/"
  counts_zip          <- paste0(local, "41598_2024_75477_MOESM2_ESM.csv.zip")
  metadata_zip        <- paste0(local, "SraRunTable (30).csv.zip")
  metadata_zip_1      <- paste0(local, "41598_2024_75477_MOESM3_ESM.csv.zip")
  metadata_zip_2      <- paste0(local, "41598_2024_75477_MOESM4_ESM.csv.zip")
  repro_counts_rds_zip<- paste0(local, "PRJNA1162476_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- paste0(local, "PRJNA1162476_dada2_taxonomy_merged.rds.zip")

>>>>>>> Stashed changes
  if (!file.exists(counts_zip)) stop("Counts file not found: ", counts_zip)

  dataset <- read.table(unz(counts_zip, "41598_2024_75477_MOESM5_ESM.csv"),
                        header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

  metadata_cols <- c("Subject", "Timepoint", "Milling", "Frequency", "StoolsperDay",
                     "BristolStoolScale_highest", "WaterContent_perc", "pH",
                     "Calprotectin_ugperg", "MPO_ngperml")

  scale_cols <- c("Subject", "Timepoint",
                  "Mean.Fungi.copies_per.mg.total.weight", 
                  "Mean.Fungi.copies_per.mg.dry.weight", 
                  "Mean.bacteria.copies_per.mg.total.weight", 
                  "Mean.bacteria.copies_per.mg.dry.weight")

  metadata <- dataset[, metadata_cols]
  scale    <- dataset[, scale_cols]
  counts_original <- dataset[, !(colnames(dataset) %in% c(metadata_cols, scale_cols))]

  taxon_names <- paste0("Taxon_", seq_len(ncol(counts_original)))
  colnames(counts_original) <- taxon_names

  tax <- tibble(
    Taxon = taxon_names,
    OriginalName = colnames(dataset)[!(colnames(dataset) %in% c(metadata_cols, scale_cols))]
  )

  counts <- bind_cols(dataset[, c("Subject", "Timepoint")], counts_original)

  row_sums <- rowSums(counts_original)
  prop_mat <- sweep(as.matrix(counts_original), 1, row_sums, FUN = "/")
  prop_mat[is.nan(prop_mat)] <- 0
  proportions <- bind_cols(dataset[, c("Subject", "Timepoint")], as_tibble(prop_mat))

  if (!file.exists(metadata_zip_1) || !file.exists(metadata_zip_2)) stop("Metadata zips not found")
  metadata1 <- read.table(unz(metadata_zip_1, "41598_2024_75477_MOESM3_ESM.csv"),
                          header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) %>%
    mutate(Replicate = gsub("Replicate ", "Replicate_", ID)) %>% select(-ID)
  metadata2 <- read.table(unz(metadata_zip_2, "41598_2024_75477_MOESM4_ESM.csv"),
                          header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) %>%
    pivot_longer(cols = starts_with("Replicate"), names_to = "Replicate", values_to = "Value") %>%
    pivot_wider(names_from = ID, values_from = Value)

  merged_metadata <- metadata1 %>%
    inner_join(metadata2, by = "Replicate") %>%
    inner_join(metadata, by = c("Subject", "Timepoint"))

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

  return(list(
    counts      = list(original = counts, reprocessed = counts_reprocessed),
    proportions = list(original = proportions, reprocessed = proportions_reprocessed),
    tax         = list(original = tax, reprocessed = tax_reprocessed),
    scale       = scale,
    metadata    = merged_metadata
  ))
}