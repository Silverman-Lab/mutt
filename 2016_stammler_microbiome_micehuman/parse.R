parse_2016_stammler_microbiome_micehuman <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
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

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2016_stammler_microbiome_micehuman")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJEB11953_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJEB11953_dada2_taxa.rds.zip")
    scale_16s_zip        <- file.path(local, "Stammler2016_scale.csv.zip")
    counts_16s_zip       <- file.path(local, "Stammler_2016_16S.csv.zip")
    metadata_16s_zip     <- file.path(local, "Stammler_2016_metadata.csv.zip")
    metadata_SRA_zip     <- file.path(local, "SraRunTable (38).csv.zip")

    # ----- Initialize everything as NA -----
    counts_original_mice <- NA
    proportions_original_mice <- NA
    tax_original_mice <- NA
    counts_original_human <- NA
    proportions_original_human <- NA
    tax_original_human <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA

    
    # ---- scale and metadata -----
    scale     <- read_zipped_table(scale_16s_zip) 
    
    # Read metadata directly with read.csv to handle # characters
    metadata_csv <- unzip(metadata_16s_zip, list = TRUE)$Name[1]
    metadata_con <- unz(metadata_16s_zip, metadata_csv)
    metadata <- read.csv(metadata_con, row.names = NULL, comment.char = "") 
    
    sra       <- read_zipped_table(metadata_SRA_zip, row.names = NULL) %>% rename(Accession = Run)
    sra <- sra %>%
      mutate(SampleID = sub(".*(MID.*)", "\\1", Sample_name))

    metadata  = merge(sra, metadata, by = "SampleID", all = TRUE)
    
    # First rename the problematic column to a simpler name
    scale <- scale %>% rename(copies = "16S rDNA copies per sample")
    
    scale = scale %>%
      left_join(metadata %>% select(SampleID, Accession), by = "SampleID") %>%
      mutate(log2_qPCR_copies = ifelse(!is.na(copies) & copies != "" & 10^copies > 0, 
                                      log2(10^copies), 
                                      NA)) %>%
      rename(log10_qPCR_copies = copies)

    # ------ original counts ------
    counts_original_mice <- read_zipped_table(counts_16s_zip)

    if (!is.na(counts_original_mice)[1]) {
      original_taxa <- colnames(counts_original_mice)

      # Create taxa mapping data frame
      tax_original_mice <- data.frame(
        Taxa = original_taxa,
        stringsAsFactors = FALSE
      )

      if (!raw) {
        align <- rename_and_align(counts_original = counts_original_mice, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
        counts_original_mice = align$counts_original
      }

      # ------ proportions from counts ------
      proportions_original_mice <- sweep(
        counts_original_mice,
        MARGIN = 1,                          
        STATS = rowSums(counts_original_mice, na.rm = TRUE),
        FUN = "/"
      )
    } else {
      proportions_original_mice <- NA
      tax_original_mice <- NA
    }

    if (file.exists(repro_counts_rds_zip)) {
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
      tax_reprocessed = make_taxa_label(tax_reprocessed)

      # ----- Convert accessions to sample IDs / Sequences to Taxa -----
      if (!raw) {
        align <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
        counts_reprocessed = align$reprocessed
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

    if (!raw) {
      counts_original_mice = fill_na_zero_numeric(counts_original_mice)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_original_mice = fill_na_zero_numeric(proportions_original_mice)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original_mice,
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = proportions_original_mice,
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = tax_original_mice,
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))
}



