parse_2017_props_isme_longitudinal <- function(raw = FALSE) {
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
  metadata_zip_2       <- file.path(local, "metadata.csv.zip")
  metadata_sra_zip     <- file.path(local, "SraRunTable (30).csv.zip")
  original_counts_zip  <- file.path(local, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list_CYANO.shared.zip")
  original_tax_zip     <- file.path(local, "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons_CYANO (1).taxonomy.zip")
  repro_counts_rds_zip <- file.path(local, "SRP066190_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "SRP066190_dada2_taxa.rds.zip")

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
  otu_matrix <- original_counts %>% select(-Sample_name) %>% as.matrix()
  otu_prop <- sweep(otu_matrix, 1, rowSums(otu_matrix), FUN = "/")
  otu_prop[is.na(otu_prop)] <- 0

  original_proportions = bind_cols(Sample = original_counts$Sample_name, as_tibble(otu_prop))

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
                dplyr::select("Sample_name", "Cell.density (cells/mL)","Cell.density.sd (cells/mL")

  meta_csv     <- unzip(metadata_sra_zip, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_sra_zip, meta_csv)
  sra          <- read.csv(meta_con) %>% as.data.frame()

  meta_csv      <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con      <- unz(metadata_zip, meta_csv)
  metadata_supp     <- read.csv(meta_con) %>% as.data.frame()
  metadata_supp     <- metadata_supp %>%
        mutate(Sample_name = as.integer(gsub("[^0-9]", "", Sample_name)))

  metadata = full_join(metadata, sra, by = "Sample_name")

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
  make_taxa_label <- function(df) {
      tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      prefixes  <- c("k", "p", "c", "o", "f", "g")
      if (!all(tax_ranks %in% colnames(df))) {
      stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
      }
      df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
      x[is.na(x) | trimws(x) == ""] <- "unclassified"
      return(x)
      })
      df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
      for (i in length(tax_ranks):1) {
          if (tax_row[i] != "unclassified") {
          return(paste0(prefixes[i], "_", tax_row[i]))
          }
      }
      return("unclassified")  
      })
      return(df)
  }
  tax_reprocessed = make_taxa_label(tax_reprocessed)

  # ----- Convert accessions to sample IDs / Sequences to Taxa -----
  # accessions to sampleIDs is study specific: IF NEED BE

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

  return(list(
    counts      = list(original = original_counts, reprocessed = counts_reprocessed),
    proportions = list(original = original_proportions, reprocessed = proportions_reprocessed),
    tax         = list(original = original_tax, reprocessed = tax_reprocessed),
    scale       = scale,
    metadata    = metadata
  ))
}