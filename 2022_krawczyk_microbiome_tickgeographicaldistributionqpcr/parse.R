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
  counts_zip           <- file.path(local, "40168_2022_1276_MOESM3_ESM.zip")
  scale_zip            <- file.path(local, "scale_qpcr.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA813158_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA813158_dada2_taxa.rds.zip")

  # ----- Scale -----
  scale_file <- unzip(scale_zip, list = TRUE)$Name[1]
  scale_con  <- unz(scale_zip, scale_file)
  scale  <- read.csv(scale_con) %>% as.data.frame() %>% mutate(rRNA_content = as.numeric(`X16S.rRNA.content.in.ng.µL`)) %>%
    mutate(log2_qPCR_16S_ng_ul = ifelse(rRNA_content > 0, log2(rRNA_content), NA)) %>%
    mutate(log10_qPCR_16S_ng_ul = ifelse(rRNA_content > 0, log10(rRNA_content), NA)) %>% 
    rename(Sample_name = `Sample.ID`) %>% select(Sample_name, rRNA_content, log2_qPCR_16S_ng_ul, log10_qPCR_16S_ng_ul)

  # ----- Metadata -----
  meta_file <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con  <- unz(metadata_zip, meta_file)
  metadata  = read.csv(meta_con) %>% as.data.frame() %>% rename(Sample_name = `Sample.ID`)
  meta_two_file <- unzip(metadata_two_zip, list = TRUE)$Name[1]
  extracted_xlsx <- unzip(metadata_two_zip, files = meta_two_file, exdir = tempdir(), overwrite = TRUE)[1]
  metadata_two <- read_excel(extracted_xlsx, sheet = 2) %>% rename(Sample_name = `Sample ID`)
  metadata <- left_join(metadata, metadata_two, by = "Sample_name") %>% 
              rename(Accession = Run)
  # -------- Counts --------
  if (file.exists(counts_zip)) {
    counts_file <- unzip(counts_zip, list = TRUE)$Name[1]
    extracted_xlsx <- unzip(counts_zip, files = counts_file, exdir = tempdir(), overwrite = TRUE)[1]
    counts_original = read_excel(extracted_xlsx, sheet = 1)
    columns_to_drop <- c("Name","Taxonomy", "Combined Abundance",	"Min",	"Max",	"Mean",	"Median",	"Std")

    # Taxa
    taxonomy_cols <- c("Sequence","Taxonomy")
    tax <- counts_original[, names(counts_original) %in% taxonomy_cols]
    tax <- tax %>%
      separate(
        col = Taxonomy,
        into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
        sep = ",\\s*",
        remove = FALSE
      ) %>%
      mutate(across(Kingdom:Genus, ~ sub("D_\\d+__", "", .)))
    tax = make_taxa_label(tax) # ----------  instead of doing this I need to reclassify the taxonomic levels with RDP for MLSCALE so for now ill just load in my other already processed file.

    # Counts
    rownames(tax) <- tax$Sequence
    counts_original <- counts_original %>% column_to_rownames("Sequence")
    counts_original = counts_original[, !(names(counts_original) %in% columns_to_drop)]
    counts_original = as.data.frame(t(counts_original))

    if (!raw) {
      aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
      counts_original = aligned$counts_original
    }

    if (!raw) {
      matched_taxa <- tax$Taxa[match(colnames(counts_original), rownames(tax))]
      colnames(counts_original) <- matched_taxa
      counts_original <- as.data.frame(t(rowsum(t(counts_original), group = colnames(counts_original))))
    }

  } else {
    counts_original = NA
  }

  proportions_original = counts_original
  proportions_original[-1] <- lapply(
    counts_original[-1],
    function(col) col / sum(col)
  )

  # ----- Reprocessed counts from RDS ZIP -----
  if (file.exists(repro_counts_rds_zip)) {
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
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
      counts_reprocessed = aligned$reprocessed
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
  }

  # DELETE LATER #####################################

  # maxwellreprocessedpreviously = file.path(local, "Krawczyk_2022_16S.csv.zip")
  # counts_original = read_zipped_table(maxwellreprocessedpreviously, row.names = NULL) %>% as.data.frame() 
  # counts_original <- counts_original %>% rownames_to_column("Sample_name") 
  # if (!raw) {
  #   aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
  #   counts_original = aligned$counts_original
  # }
  # rownames(counts_original) <- counts_original$Sample_name
  # counts_original <- counts_original[, -1, drop = FALSE]  
  # proportions_original <- sweep(counts_original, MARGIN = 1,STATS  = rowSums(counts_original), FUN = "/")
  ####################################################
  if (!raw) {
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

  # ----- Return all -----
  return(list(
    counts      = list(
      original = counts_original,
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = proportions_reprocessed
    ),
    tax         = list(
      original = tax,
      reprocessed = tax_reprocessed
    ),
    scale       = scale,
    metadata    = metadata
  ))
}
