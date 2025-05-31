<<<<<<< Updated upstream
parse_2020_tettamantiboshier_msystems_vaginaltimeseries <- function(raw = FALSE, align = FALSE) {
=======
parse_2020_tettamantiboshier_msystems_vaginaltimeseries <- function() {
>>>>>>> Stashed changes
  required_pkgs <- c("tibble", "tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
<<<<<<< Updated upstream
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

  counts_original <- NA
  proportions_original <- NA
  tax_original <- NA
  scale <- NA
  metadata <- NA
  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA
  counts_reprocessed2 <- NA
  proportions_reprocessed2 <- NA
  tax_reprocessed2 <- NA

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

    metadata = metadata %>% full_join(scale %>% select(Sample_ID, Participant, Hours_In_Study, `Lactobacillus_crispatus`, `Lactobacillus_jensenii`, `Lactobacillus_iners`, `Gardnerella_vaginalis`, `Megasphaera`, `BVAB2`, `Atopobium_vaginae`, `Accession`), by = "Accession")
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

    proportions_original <- sweep(counts_original, 1, rowSums(counts_original), FUN = "/")

  } else {
    counts_original      <- NA
    proportions_original <- NA
    tax_original         <- NA
  }

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
    temp_dir <- tempdir("repro")
    dir.create(temp_dir)
    unzipped = unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
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
    unzipped = unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
        aligned_data = rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample_ID", align = align, study_name = basename(local))
        counts_reprocessed = aligned_data$reprocessed
        counts_reprocessed2 = aligned_data$reprocessed
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
    proportions_reprocessed <- as.data.frame(proportions_reprocessed, check.names = FALSE)
    cleanup_tempfiles(temp_dir)
  }

  if (!raw) {
    counts_original = fill_na_zero_numeric(counts_original)
    counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
    proportions_original = fill_na_zero_numeric(proportions_original)
    proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
  }
=======
  library(tibble)
  library(tidyverse)

  # ----- File paths -----
  local               <- file.path("/2020_tettamantiboshier_msystems_vaginaltimeseries/")
  counts_zip          <- paste0(local, "originalcounts.csv.zip")
  scale_zip           <- paste0(local, "scale.csv.zip")
  metadata_zip        <- paste0(local, "SraRunTable (30).csv.zip")
  repro_counts_rds_zip<- paste0(local, "PRJNA549339_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- paste0(local, "PRJNA549339_dada2_taxonomy_merged.rds.zip")

  # ----- Read counts & proportions -----
  if (file.exists(counts_zip)) {
    counts_file <- unzip(counts_zip, list = TRUE)$Name[1]
    counts_con  <- unz(counts_zip, counts_file)
    countsdata  <- read.csv(counts_con) %>% as.data.frame()

    columns_to_drop <- c("Participant", "Hours_In_Study")
    counts_original <- countsdata[, !(names(countsdata) %in% columns_to_drop)]
    metadatacols    <- countsdata[, names(countsdata) %in% columns_to_drop]

    taxon_names <- paste0("Taxon_", seq_len(ncol(counts_original)))
    colnames(counts_original) <- taxon_names

    tax <- tibble(
      Taxon = taxon_names,
      OriginalName = colnames(countsdata)[!(colnames(countsdata) %in% columns_to_drop)]
    )

    counts <- bind_cols(metadatacols, counts_original)

    row_sums    <- rowSums(counts_original)
    prop_mat    <- sweep(as.matrix(counts_original), 1, row_sums, FUN = "/")
    prop_mat[is.nan(prop_mat)] <- 0
    proportions <- bind_cols(metadatacols, as_tibble(prop_mat))

  } else {
    counts      <- NA
    proportions <- NA
    tax         <- NA
  }

  # ----- Scale data -----
  if (file.exists(scale_zip)) {
    scale_file <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_con  <- unz(scale_zip, scale_file)
    scale  = read.csv(scale_con) %>% as.data.frame()

  } else {
    scale <- NA
  }

  # ----- Metadata -----
  if (file.exists(metadata_zip)) {
    metadata_file <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_con  <- unz(metadata_zip, metadata_file)
    metadata = read.csv(metadata_con) %>% as.data.frame()
    # This needs to be cleaned up because im not sure how to link the sample IDs.
  } else {
    metadata <- NA
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
>>>>>>> Stashed changes

  # ----- Return combined object -----
  return(list(
    counts      = list(
<<<<<<< Updated upstream
      original = counts_original,
      reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
      ),
    proportions = list(
      original = proportions_original,
      reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
      ),
    tax         = list(
      original = tax_original,
      reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
=======
      original = counts,
      reprocessed = counts_reprocessed
      ),
    proportions = list(
      original = proportions,
      reprocessed = proportions_reprocessed
      ),
    tax         = list(
      original = tax,
      reprocessed = tax_reprocessed
>>>>>>> Stashed changes
      ),
    scale       = scale,
    metadata    = metadata
  ))
}
