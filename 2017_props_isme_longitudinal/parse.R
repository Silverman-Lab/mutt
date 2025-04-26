parse_2017_props_isme_longitudinal <- function() {
  required_pkgs <- c("tibble", "tidyverse", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)
  library(readr)

    # ----- Local base directory -----
    local <- file.path("2017_props_isme_longitudinal")

    # ----- File paths -----
    metadata_zip         <- file.path(local, "supplemental_metadata.csv.zip")
    metadata_sra_zip     <- file.path(local, "SraRunTable (30).csv.zip")
    original_counts_zip  <- file.path(local, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list_CYANO.shared.zip")
    original_tax_zip     <- file.path(local, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons_CYANO (1).taxonomy.zip")
    repro_counts_rds_zip <- file.path(local, "SRP066190_dada2_merged_nochim.rds.zip")
    repro_tax_zip        <- file.path(local, "SRP066190_dada2_taxonomy_merged.rds.zip")

  # ------ original counts, proportions, tax ---- 

    shared_file <- unzip(original_counts_zip, list = TRUE)$Name[1]
    shared_con  <- unz(original_counts_zip, shared_file)
    shared_raw <- read.table(shared_con, header = TRUE, sep = "\t", check.names = FALSE)

    otu_cols <- names(shared_raw)[-(1:3)]
    taxon_names <- paste0("Taxon_", seq_along(otu_cols))
    taxonomy_key <- tibble(Taxon = taxon_names, OriginalOTU = otu_cols)

  original_counts = shared_raw %>%
    select(-label, -numOtus) %>%
    rename(Sample_name = Group) %>%
    setNames(c("Sample_name", taxon_names))

    tax_file <- unzip(original_tax_zip, list = TRUE)$Name[1]
    tax_con  <- unz(original_tax_zip, tax_file)
    taxonomy_df <- read.table(tax_con, sep = "\t", header = FALSE,
                            col.names = c("OriginalOTU", "Taxonomy", "Confidence"),
                            stringsAsFactors = FALSE)

    taxonomy_mapped <- left_join(taxonomy_key, taxonomy_df, by = "OriginalOTU")

  original_tax = taxonomy_mapped %>%
    separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order",
                                "Family", "Genus", "Species"),
            sep = ";", fill = "right", extra = "drop") %>%
    mutate(across(Kingdom:Species, str_trim)) %>%
    select(Taxon, Kingdom, Phylum, Class, Order, Family, Genus, Species, Confidence)

    otu_matrix <- original_counts %>% select(-Sample) %>% as.matrix()
    otu_prop <- sweep(otu_matrix, 1, rowSums(otu_matrix), FUN = "/")
    otu_prop[is.na(otu_prop)] <- 0

  original_proportions = bind_cols(Sample = original_counts$Sample, as_tibble(otu_prop))

  # ------ metadata and scale -----
  meta_csv     <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_zip, meta_csv)
  metadata     <- read.csv(meta_con) %>% as.data.frame()
  metadata     <- metadata %>%
        mutate(Sample_name = as.integer(gsub("[^0-9]", "", Sample_name)))

  scale = metadata %>%
                dplyr::select("Sample_name", "Cell.density (cells/mL)","Cell.density.sd (cells/mL")

  meta_csv     <- unzip(metadata_sra_zip, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_sra_zip, meta_csv)
  sra          <- read.csv(meta_con) %>% as.data.frame()

  metadata = full_join(metadata, sra, by = "Sample_name")

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

  sample_map <- metadata %>%
    filter(Run %in% colnames(counts_reprocessed)) %>%
    distinct(Run, Sample_name)

  relabel <- function(df, map_df) {
    rename_vec <- setNames(map_df$Sample_name, map_df$Run)
    colnames(df)[colnames(df) %in% names(rename_vec)] <- rename_vec[colnames(df)[colnames(df) %in% names(rename_vec)]]
    df
  }

  counts_reprocessed      <- relabel(counts_reprocessed, sample_map)
  proportions_reprocessed <- relabel(proportions_reprocessed, sample_map)

  return(list(
    counts      = list(original = original_counts, reprocessed = counts_reprocessed),
    proportions = list(original = original_proportions, reprocessed = proportions_reprocessed),
    tax         = list(original = original_tax, reprocessed = tax_reprocessed),
    scale       = scale,
    metadata    = metadata
  ))
}