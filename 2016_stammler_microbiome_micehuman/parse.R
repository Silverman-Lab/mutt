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
    counts_reprocessed2 <- NA
    proportions_reprocessed2 <- NA
    tax_reprocessed2 <- NA

    
    # ---- scale and metadata -----
    scale     <- read_zipped_table(scale_16s_zip) 
    
    # Read metadata directly with read.csv to handle # characters
    metadata_csv <- unzip(metadata_16s_zip, list = TRUE)$Name[1]
    metadata_con <- unz(metadata_16s_zip, metadata_csv)
    metadata <- read.csv(metadata_con, row.names = NULL, comment.char = "") 
    
    sra       <- read_zipped_table(metadata_SRA_zip, row.names = NULL) %>% rename(Accession = Run)
    sra <- sra %>%
      mutate(SampleID = sub(".*(MID.*)", "\\1", Sample_name))

    metadata  = full_join(sra, metadata, by = "SampleID")
    
    # First rename the problematic column to a simpler name
    scale <- scale %>% rename(copies = "16S rDNA copies per sample")
    
    scale = scale %>%
      left_join(metadata %>% select(SampleID, Accession), by = "SampleID") %>%
      mutate(log2_qPCR_copies = ifelse(!is.na(copies) & copies != "" & copies != 0, copies * log2(10), NA)) %>%
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
        aligned <- rename_and_align(counts_original = counts_original_mice, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
        counts_original_mice = aligned$counts_original
        original_names <- colnames(counts_original_mice)
        counts_original_mice <- as.data.frame(lapply(counts_original_mice, as.numeric), row.names = rownames(counts_original_mice), col.names = original_names, check.names = FALSE)
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

      # ----- Taxonomy reprocessed ----
      unzipped = unzip(repro_tax_zip, exdir = temp_rds, overwrite = TRUE)
      tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
      if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
      tax_reprocessed <- as.data.frame(readRDS(tax_file))

      # ----- Convert sequences to lowest rank taxonomy found and update key -----
      tax_reprocessed = make_taxa_label(tax_reprocessed)

      # ----- Convert accessions to sample IDs / Sequences to Taxa -----
      if (!raw) {
        aligned <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
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
      counts_original_mice = fill_na_zero_numeric(counts_original_mice)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_original_mice = fill_na_zero_numeric(proportions_original_mice)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
      proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    }

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original_mice,
            reprocessed = list(
              rdp19 = counts_reprocessed,
              rdp16 = counts_reprocessed2
            )
        ),
        proportions = list(
            original = proportions_original_mice,
            reprocessed = list(
              rdp19 = proportions_reprocessed,
              rdp16 = proportions_reprocessed2
            )
        ),
        tax = list(
            original = tax_original_mice,
            reprocessed = list(
              rdp19 = tax_reprocessed,
              rdp16 = tax_reprocessed2
            )
        ),
        scale = scale,
        metadata = metadata
    ))
}


