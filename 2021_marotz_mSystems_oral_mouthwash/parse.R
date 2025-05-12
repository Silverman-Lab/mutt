parse_2021_marotz_mSystems_oral_mouthwash <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  if (!is.logical(raw)) {
    stop("raw must be a logical value")
  }
  if (!is.logical(align)) {
    stop("align must be a logical value")
  }

  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2021_marotz_mSystems_oral_mouthwash")

  # ----- File paths -----
  counts_zip    <- file.path(local, "2021_marotz_mSystems_oral_mouthwash.RDS.zip")
  metadata_zip  <- file.path(local, "T3_SRS_metadata_ms.txt.zip")
  sra_zip       <- file.path(local, "SraRunTable (40).csv.zip")

  repro_counts_zips <- c(
    file.path(local, "ERP111447_dada2_counts.rds.zip"),
    file.path(local, "ERP117149_dada2_counts.rds.zip")
  )

  repro_tax_zips <- c(
    file.path(local, "ERP111447_dada2_taxa.rds.zip"),
    file.path(local, "ERP117149_dada2_taxa.rds.zip")
  )

  # ----- Metadata and Scale -----
  metadata_txt <- unzip(metadata_zip, list = TRUE)$Name[1]
  metadata <- read.csv(unz(metadata_zip, metadata_txt), sep = ",")
  sra <- read_zipped_table(sra_zip, row.names = NULL) %>% rename(Accession = Run, saliva_sample_ID = saliva_sample_id) 
  sra$SampleID <- paste0(sra$saliva_sample_ID, ".", sra$timepoint, ".", sra$processing)
  metadata$Sample <- paste0(metadata$saliva_sample_ID, ".", metadata$processing)
  metadata <- merge(sra, metadata, by = "SampleID")
  
  # Add replicate column based on run_date within SampleID groups
  metadata <- metadata %>%
    arrange(Sample, `run_date (exp)`) %>%
    group_by(Sample) %>%
    mutate(replicate = row_number()) %>%
    ungroup()

  metadata$SampleID <- paste0(metadata$Sample, ".", metadata$replicate)
    
  metadata = remove_empty_columns(metadata)

  scale <- metadata %>% select(SampleID, FC_cells_per_ul_r1, FC_cells_per_ul_r2, FC_avg_cells_per_ul, 
                              FC_avg_cells_5_min, qpcr_median_16s_copies_per_2ul_dna) %>% 
           mutate(FC_sd_cells_per_ul = sqrt((FC_cells_per_ul_r1 - FC_avg_cells_per_ul)^2 + (FC_cells_per_ul_r2 - FC_avg_cells_per_ul)^2) / 2)

  scale = scale %>% 
    mutate(log2_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log2(FC_avg_cells_per_ul), NA)) %>%
    mutate(log10_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log10(FC_avg_cells_per_ul), NA)) %>%
    mutate(log2_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log2(FC_avg_cells_5_min), NA)) %>%
    mutate(log10_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log10(FC_avg_cells_5_min), NA)) %>%
    mutate(log2_FC_sd_cells_per_ul = ifelse(FC_sd_cells_per_ul > 0, log2(FC_sd_cells_per_ul), NA)) %>%
    mutate(log10_FC_sd_cells_per_ul = ifelse(FC_sd_cells_per_ul > 0, log10(FC_sd_cells_per_ul), NA)) %>%
    mutate(log2_qpcr_median_16s_copies_per_2ul_dna = ifelse(qpcr_median_16s_copies_per_2ul_dna > 0, log2(qpcr_median_16s_copies_per_2ul_dna), NA)) %>%
    mutate(log10_qpcr_median_16s_copies_per_2ul_dna = ifelse(qpcr_median_16s_copies_per_2ul_dna > 0, log10(qpcr_median_16s_copies_per_2ul_dna), NA))

  # ----- Original Counts and Taxonomy -----
  orig_rds_file <- unzip(counts_zip, list = TRUE)$Name[1]
  temp_dir <- tempfile("repro")
  dir.create(temp_dir)
  unzip(counts_zip, files = orig_rds_file, exdir = temp_dir, overwrite = TRUE)
  counts_original <- readRDS(file.path(temp_dir, orig_rds_file)) %>% t() %>% as.data.frame()
  counts_original <- collapse_duplicate_columns_exact(counts_original)
  original_names <- colnames(counts_original)
  counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
  cleanup_tempfiles(temp_dir)
  if (!raw) {
    aligned = rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
    counts_original = aligned$counts_original
    original_names <- colnames(counts_original)
    counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
  }
  proportions_original <- sweep(counts_original, 1, rowSums(counts_original), "/")


  raw_tax <- data.frame(Taxa = colnames(counts_original))
  tax_original <- raw_tax %>%
    mutate(taxa = str_trim(Taxa)) %>%
    separate(
      Taxa,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = "\\s*;\\s*",
      extra = "drop",
      fill = "right"
    )

  tax_original = make_taxa_label(tax_original)

  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA

  if (all(file.exists(repro_counts_zips), file.exists(repro_tax_zips))) {
    # # Process multiple zipped RDS files
    counts_reprocessed_list <- list()
    proportions_reprocessed_list <- list()
    tax_reprocessed_list <- list()

    for (i in seq_along(repro_counts_zips)) {
        # Unzip and read counts
        temp_rds <- tempfile("repro")
        dir.create(temp_rds)
        unzipped = unzip(repro_counts_zips[i], exdir = temp_rds, overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop(paste("No *_counts.rds file found for index", i))
        counts <- as.data.frame(readRDS(counts_file))

        # Unzip and read taxonomy
        unzipped = unzip(repro_tax_zips[i], exdir = temp_rds, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop(paste("No *_taxa.rds file found for index", i))
        tax <- as.data.frame(readRDS(tax_file))
        tax <- make_taxa_label(tax)

        if (!raw) {
          aligned = rename_and_align(counts_reprocessed = counts, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
          counts = aligned$reprocessed
          matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
          colnames(counts) <- matched_taxa
          counts = collapse_duplicate_columns_exact(counts)
          original_names <- colnames(counts)
          counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
        }
        # proportions
        proportions <- sweep(counts, 1, rowSums(counts), '/')

        # Label with study name based on zip filename prefix
        study_id <- sub("_.*$", "", basename(tools::file_path_sans_ext(repro_counts_zips[i])))
        counts$Study <- study_id
        proportions$Study <- study_id
        tax$Study <- study_id

        counts_reprocessed_list[[i]] <- counts
        proportions_reprocessed_list[[i]] <- proportions
        tax_reprocessed_list[[i]] <- tax

        cleanup_tempfiles(temp_rds)
    }

    # Combine all
    counts_reprocessed <- bind_rows(counts_reprocessed_list)
    proportions_reprocessed <- bind_rows(proportions_reprocessed_list)
    tax_reprocessed <- bind_rows(tax_reprocessed_list)

  } else {
  # DELETE LATER WHEN REPROCESS HAS FINALIZED:
  counts_zip = file.path(local, "Marotz_2021_16S.csv.zip")
  counts_reprocessed = read_zipped_table(counts_zip) %>% as.data.frame()
  counts_reprocessed <- counts_reprocessed %>% mutate(across(everything(), as.numeric))
  if (!raw) {
    aligned = rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
    counts_reprocessed = aligned$counts_original
    original_names <- colnames(counts_reprocessed)
    counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
  }
  proportions_reprocessed <- sweep(counts_reprocessed, MARGIN = 1,STATS  = rowSums(counts_reprocessed), FUN = "/")
  tax_reprocessed <- data.frame(Taxa = colnames(counts_reprocessed), stringsAsFactors = FALSE)
  ############################################
  }

  if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      proportions_original = fill_na_zero_numeric(proportions_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

  # ----- Return structured list -----
  return(list(
    counts = list(
      original = counts_original,
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = proportions_reprocessed
    ),
    tax = list(
      original = tax_original,
      reprocessed = tax_reprocessed
    ),
    scale = scale,
    metadata = metadata
  ))
}
