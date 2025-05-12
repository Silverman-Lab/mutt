parse_2020_tettamantiboshier_msystems_vaginaltimeseries <- function(raw = FALSE, align = FALSE) {
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
  local <- file.path("2020_tettamantiboshier_msystems_vaginaltimeseries")

  # ----- File paths -----
  counts_zip           <- file.path(local, "originalcounts.csv.zip")
  scale_zip            <- file.path(local, "scale_connect20250425.csv.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA549339_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA549339_dada2_taxa.rds.zip")

  # ----- Scale data -----
  if (file.exists(scale_zip)) {
    scale_file <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_con  <- unz(scale_zip, scale_file)
    scale  = read.csv(scale_con) %>% as.data.frame() %>% rename(Accession = Run, Sample = `Sample.Name`, Sample_ID = `X`) 
    scale = scale %>% 
      mutate(log2_total_16S = log2(Total_16S),
             log10_total_16S = log10(Total_16S)) %>%
      mutate(across(c(log2_total_16S, log10_total_16S), ~ifelse(is.infinite(.), NA, .)))
  } else {
    scale <- NA
  }

  # ----- Metadata -----
  if (file.exists(metadata_zip)) {
    metadata_file <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_con  <- unz(metadata_zip, metadata_file)
    metadata = read.csv(metadata_con) %>% as.data.frame() %>% rename(Accession = Run, Sample = `Sample.Name`)
  } else {
    metadata <- NA
  }

   # ----- Read counts & proportions -----
  if (file.exists(counts_zip)) {
    countsdata  <- read_zipped_table(counts_zip) 

    counts_original <- countsdata[, !(names(countsdata) %in% c("Participant", "Hours_In_Study"))]
    metadatacols    <- countsdata[, c("Participant", "Hours_In_Study") ] %>% rownames_to_column(var = "Sample_ID")

    metadata = metadata %>% merge(scale %>% select(Sample_ID, Participant, Hours_In_Study, `Lactobacillus_crispatus`, `Lactobacillus_jensenii`, `Lactobacillus_iners`, `Gardnerella_vaginalis`, `Megasphaera`, `BVAB2`, `Atopobium_vaginae`, `Accession`), by = "Accession")
    scale = scale %>% select(Sample_ID, Participant, Hours_In_Study, Accession, `log2_total_16S`, `log10_total_16S`)

    tax_original <- tibble(
      Taxa = colnames(countsdata)[!(colnames(countsdata) %in% c("Participant", "Hours_In_Study"))]
    )

    if (!raw) {
      aligned_data = rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "Sample_ID", align = align, study_name = basename(local))
      counts_original = aligned_data$counts_original
      original_names <- colnames(counts_original)
      counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
    }

    prop_mat    <- sweep(as.matrix(counts_original), 1, rowSums(counts_original), FUN = "/")
    proportions_original <- bind_cols(metadatacols, as_tibble(prop_mat))

  } else {
    counts_original      <- NA
    proportions_original <- NA
    tax_original         <- NA
  }

  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA 

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
    temp_dir <- tempdir("repro")
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
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
        aligned_data = rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample_ID", align = align, study_name = basename(local))
        counts_reprocessed = aligned_data$reprocessed
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
        original_names <- colnames(counts_reprocessed)
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
    }
    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
    proportions_reprocessed <- as.data.frame(proportions_reprocessed, check.names = FALSE)
    cleanup_tempfiles(temp_dir)
  }


  if (!raw) {
    counts_original = fill_na_zero_numeric(counts_original)
    counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
    proportions_original = fill_na_zero_numeric(proportions_original)
    proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

  # ----- Return combined object -----
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
      original = tax_original,
      reprocessed = tax_reprocessed
      ),
    scale       = scale,
    metadata    = metadata
  ))
}
