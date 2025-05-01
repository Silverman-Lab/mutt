parse_2022_krawczyk_microbiome_tickgeographicaldistributionqpcr <- function() {
  required_pkgs <- c("tibble", "tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)
  library(readxl)

<<<<<<< Updated upstream
  # ----- Local base directory -----
  local <- file.path("2022_krawczyk_microbiome_tickgeographicaldistributionqpcr")

  # ----- File paths -----
  metadata_two_zip     <- file.path(local, "40168_2022_1276_MOESM1_ESM.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv")
  counts_zip           <- file.path(local, "40168_2022_1276_MOESM3_ESM.zip")
  scale_zip            <- file.path(local, "scale_qpcr.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA813158_dada2_merged_nochim.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA813158_dada2_taxonomy_merged.rds.zip")

=======
  # ----- File paths -----
  local               <- "2022_krawczyk_microbiome_tickgeographicaldistributionqpcr/"
  metadata_two_zip    <- paste0(local, "40168_2022_1276_MOESM1_ESM.zip")
  metadata_zip        <- paste0(local, "SraRunTable (30).csv")
  counts_zip          <- paste0(local, "40168_2022_1276_MOESM3_ESM.zip")
  scale_zip           <- paste0(local, "scale_qpcr.zip")
  repro_counts_rds_zip<-  paste0(local, "PRJNA813158_dada2_merged_nochim.rds.zip"
  repro_tax_zip       <-  paste0(local, "PRJNA813158_dada2_taxonomy_merged.rds.zip"

>>>>>>> Stashed changes
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
    scale       = scale,
    metadata    = metadata
  ))
}
