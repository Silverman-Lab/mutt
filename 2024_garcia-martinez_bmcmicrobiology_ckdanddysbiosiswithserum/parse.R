parse_2024_garciamartinez_bmcmicrobiology_ckdflow <- function() {
    
    required_pkgs <- c("tidyverse", "readxl")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    
    if (length(missing_pkgs) > 0) {
      stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
           ". Please install them before running this function.")
    }
    
    # Load needed libraries 
    library(tidyverse)
    library(readxl)
    
    ## Define ZIP file paths
    metadata_zip <- "2024_garcia-martinez_bmcmicrobiology_ckdanddysbiosiswithserum/SraRunTable (19).csv.zip"
    scale_zip <- "2024_garcia-martinez_bmcmicrobiology_ckdanddysbiosiswithserum/12866_2024_3590_MOESM1_ESM.xlsx.zip"
    orig_counts_zip     <- "2017_vandeputte_nature_flow/OTU_nochim.zip"
    orig_tax_rdp_zip    <- NA
    orig_tax_silva_zip  <- NA
    orig_prop_zip       <- NA
    repro_counts_rds_zip<- "2024_garcia-martinez_bmcmicrobiology_ckdanddysbiosiswithserum/PRJEB67373_dada2_merged_nochim.rds.zip"
    repro_tax_zip       <- "2024_garcia-martinez_bmcmicrobiology_ckdanddysbiosiswithserum/PRJEB67373_dada2_taxonomy_merged.rds.zip"

    ## Read Metadata
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_con <- unz(metadata_zip, metadata_csv)
    metadata <- read.csv(metadata_con, row.names = "Sample.Name") %>%
      as.data.frame()
    
    ## Read Scale
    scale_files <- unzip(scale_zip, list = TRUE)
    scale_xlsx <- scale_files$Name[1]
    temp_dir <- tempdir()
    unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
    scale_path <- file.path(temp_dir, scale_xlsx)
    
    scale_raw <- read_xlsx(scale_path, sheet = 1) %>%
      separate(col = 1, into = c("population", "Sample.Name"), sep = "_", remove = TRUE)
    metadata <- metadata %>%
      rownames_to_column("Sample.Name") %>%
      left_join(scale_raw %>% select(Sample.Name, population), by = "Sample.Name") %>%
      column_to_rownames("Sample.Name")
    scale <- scale_raw %>%
      select(-population) %>%
      column_to_rownames("Sample.Name")

    ## Original Counts
    if (file.exists(orig_counts_zip)) {
      orig_csv <- unzip(orig_counts_zip, list = TRUE)$Name[1]
      orig_con <- unz(orig_counts_zip, orig_csv)
      orig_mat <- read.csv(orig_con, row.names = 1, check.names = FALSE)
      counts_original <- as.data.frame(orig_mat)
      counts_original$Sequence <- rownames(counts_original)
      counts_original <- counts_original[, c("Sequence", setdiff(names(counts_original), "Sequence"))]
      rownames(counts_original) <- paste0("Taxon_", seq_len(nrow(counts_original)))
    } else {
      counts_original <- NA
    }

    # Original Proportions
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

    # Original taxonomies 
    read_taxonomy_zip <- function(zip_path) {
      if (!file.exists(zip_path)) return(NA)
      zip_contents <- unzip(zip_path, list = TRUE)
      csv_name <- zip_contents$Name[1]
      con <- unz(zip_path, csv_name)
      tax_df <- read.csv(con, row.names = 1, check.names = FALSE)
      tax_df$Sequence <- rownames(tax_df)
      tax_df <- tax_df[, c("Sequence", setdiff(names(tax_df), "Sequence"))]
      rownames(tax_df) <- paste0("Taxon_", seq_len(nrow(tax_df)))
      as_tibble(tax_df, rownames = "Taxon")
    }

    tax_original_rdp   <- read_taxonomy_zip(orig_tax_rdp_zip)
    tax_original_silva <- read_taxonomy_zip(orig_tax_silva_zip)

    ## Taxonomy reprocessed 
    temp_tax <- tempfile(fileext = ".rds")
    unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
    tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
    taxonomy_matrix <- readRDS(tax_file)
    rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
    tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
    tax_reprocessed <- tax_table


    return(list(
      counts      = list(
        original    = counts_original,
        reprocessed = counts_reprocessed
      ),
      proportions = list(
        original    = proportions_original,
        reprocessed = proportions_reprocessed
      ),
      tax         = list(
        original_rdp    = tax_original_rdp,
        original_silva  = tax_original_silva,
        reprocessed = tax_reprocessed
      ),
      scale       = scale,
      metadata    = metadata
    ))
}
  