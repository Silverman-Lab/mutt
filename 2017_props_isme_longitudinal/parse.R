parse_2017_props_isme_longitudinal <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "readr")
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
  library(tibble)
  library(tidyverse)
  library(readr)

  # ----- Local base directory -----
  local <- file.path("2017_props_isme_longitudinal")

  # ----- File paths -----
  metadata_zip         <- file.path(local, "supplemental_metadata.csv.zip")
  metadata_zip_2       <- file.path(local, "metadata.csv.zip")
  metadata_sra_zip     <- file.path(local, "SraRunTable (30).csv.zip")
  original_counts_zip  <- file.path(local, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list_CYANO.shared.zip")
  original_tax_zip     <- file.path(local, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons_CYANO (1).taxonomy.zip")
  repro_counts_rds_zip <- file.path(local, "SRP066190_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "SRP066190_dada2_taxa.rds.zip")

  # ------ metadata and scale -----
  meta_csv     <- unzip(metadata_zip_2, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_zip_2, meta_csv)
  metadata     <- read.csv(meta_con) %>% as.data.frame()
  metadata     <- metadata %>%
        mutate(Sample_name = as.integer(gsub("[^0-9]", "", Sample_name)))
  metadata <- metadata %>%
    dplyr::rename(
      `Cell.density (cells/mL)` = Cell.density..cells.mL.,
      `Cell.density.sd (cells/mL)` = Cell.density.sd..cells.mL.
    )
  scale = metadata %>%
                dplyr::select("Sample_name", "Cell.density (cells/mL)","Cell.density.sd (cells/mL") %>%
                mutate(log2_FC_cell_ml = ifelse(`Cell.density (cells/mL)`>0, log2(`Cell.density (cells/mL)`), NA)) %>%
                mutate(log10_FC_cell_ml = ifelse(`Cell.density (cells/mL)`>0, log10(`Cell.density (cells/mL)`), NA))

  meta_csv     <- unzip(metadata_sra_zip, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_sra_zip, meta_csv)
  sra          <- read.csv(meta_con) %>% as.data.frame()

  meta_csv      <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con      <- unz(metadata_zip, meta_csv)
  metadata_supp     <- read.csv(meta_con) %>% as.data.frame()
  metadata_supp     <- metadata_supp %>%
        mutate(Sample_name = as.integer(gsub("[^0-9]", "", Sample_name)))

  metadata = full_join(metadata, sra, by = "Sample_name")

  # ------ original counts, proportions, tax ---- 
  shared_file <- unzip(original_counts_zip, list = TRUE)$Name[1]
  shared_con  <- unz(original_counts_zip, shared_file)
  shared_raw <- read.table(shared_con, header = TRUE, sep = "\t", check.names = FALSE)

  original_counts = shared_raw %>%
    select(-label, -numOtus) %>%
    rename(Sample_name = Group) 

  tax_file <- unzip(original_tax_zip, list = TRUE)$Name[1]
  tax_con  <- unz(original_tax_zip, tax_file)
  taxonomy <- read.table(
    tax_con,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    fill = TRUE,
    row.names = NULL,   
    check.names = FALSE    
  )
  if (!raw) {
    aligned <- rename_and_align(counts_original = original_counts, metadata = metadata, scale = scale, by_col = "Sample_name", align = align, study_name = basename(local))
    original_counts = aligned$counts_original
  }
  otu_matrix <- original_counts %>% select(-Sample_name)
  otu_prop <- sweep(otu_matrix, 1, rowSums(otu_matrix), FUN = "/")

  original_proportions = bind_cols(Sample = original_counts$Sample_name, as_tibble(otu_prop))

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
    temp_rds <- tempfile("repro")
    dir.create(temp_rds)
    unzipped = unzip(repro_counts_rds_zip, exdir = temp_rds, overwrite = TRUE)
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
        NULL
      }
       
      } else {
        tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
    }

    # ----- Taxonomy reprocessed -----
    unzipped = unzip(repro_tax_zip, exdir = temp_rds, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))

    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample_name", align = align, study_name = basename(local))
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
    cleanup_tempfiles(temp_rds)
  }
  
  if (!raw) {
    original_counts = fill_na_zero_numeric(original_counts)
    counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
    original_proportions = fill_na_zero_numeric(original_proportions)
    proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
    proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
  }

  return(list(
    counts      = list(original = original_counts, reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)),
    proportions = list(original = original_proportions, reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)),
    tax         = list(original = original_tax, reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)),
    scale       = scale,
    metadata    = metadata
  ))
}