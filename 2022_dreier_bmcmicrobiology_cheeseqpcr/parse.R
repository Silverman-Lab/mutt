parse_2022_dreier_bmcmicrobiology_cheeseqpcr <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "readxl")
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
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2022_dreier_bmcmicrobiology_cheeseqpcr")

  # ----- File paths -----
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  scale_zip            <- file.path(local, "V18-22-21_htqpcr_rawdata.csv.zip")
  orig_counts_zip      <- file.path(local, "V18-22-21_ASV_counts_table.csv.zip")
  orig_prop_zip        <- NA
  repro_counts_rds_zip <- file.path(local, "PRJNA786903_dada2_counts.rds.zip") # Taxa need to be identified using cheese database
  repro_tax_zip        <- file.path(local, "PRJNA786903_dada2_taxa.rds.zip") # Taxa need to be identified using cheese database
  deletewhenfinished   <- file.path(local, "Dreier_2022_16S.csv.zip")

  # Initialize variables
  counts_original <- NA
  proportions_original <- NA
  tax_original <- NA
  counts_original2 <- NA
  proportions_original2 <- NA
  tax_original2 <- NA
  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA
  counts_reprocessed2 <- NA
  proportions_reprocessed2 <- NA
  tax_reprocessed2 <- NA
  scale <- NA
  metadata <- NA

  # Scale
  if (file.exists(scale_zip)) {
    temp_dir <- tempfile("scale")
    dir.create(temp_dir)
    # Extract raw Fluidigm CSV from the zip archive
    scale_file <- unzip(scale_zip, list = TRUE)$Name[1]
    csv_file   <- unzip(scale_zip, files = scale_file, exdir = temp_dir, overwrite = TRUE)[1]

    # Read and flatten Fluidigm 2-row header
    raw_hdr <- read_lines(csv_file, n_max = 12)
    hdr1    <- str_split(raw_hdr[11], ",")[[1]] %>% str_trim()
    hdr2    <- str_split(raw_hdr[12], ",")[[1]] %>% str_trim()
    col_nms <- make.names(paste0(hdr1, "_", hdr2), allow_ = TRUE)

    # Load qPCR data
    df <- read_csv(csv_file, skip = 12, col_names = col_nms, show_col_types = FALSE) 
    df <- df %>%
      mutate(
        Ct_Value            = na_if(as.character(Ct_Value), "Undetermined"),
        Ct_Value            = as.numeric(Ct_Value),
        Sample_rConc        = as.numeric(Sample_rConc),
        Ct_Calibrated.rConc = as.numeric(Ct_Calibrated.rConc)
      )
    # 1. Total abundance from calibrated concentrations
    df %>%
        group_by(Sample_Name) %>%
        summarise(
        total_copies = sum(Ct_Calibrated.rConc, na.rm = TRUE),
        mean_copies  = mean(Ct_Calibrated.rConc, na.rm = TRUE),
        sd_copies    = sd(Ct_Calibrated.rConc, na.rm = TRUE),
        n_reps       = sum(!is.na(Ct_Calibrated.rConc)),
        .groups = "drop"
        ) 
    # 2. Fit standard curve from "Standard" wells
    stds <- df %>%
        filter(
        str_to_lower(Sample_Type) == "standard",
        !is.na(Ct_Value),
        !is.na(Sample_rConc)
        ) %>%
        mutate(log_copies = log10(Sample_rConc))

    if (nrow(stds) < 2) {
        stop("✗ Not enough standard wells found to fit a Ct → copy curve.")
    }

    fit <- lm(Ct_Value ~ log_copies, data = stds)
    m   <- coef(fit)["log_copies"]
    b   <- coef(fit)["(Intercept)"]
    #message(sprintf("✓ Standard curve: Ct = %.4f × log10(copies) + %.4f", m, b))

    # 3. Apply standard curve to unknowns
    scale <- df %>%
      filter(
        !is.na(Ct_Value),
        str_to_lower(Sample_Type) != "standard" & str_to_lower(Sample_Name) != "not relevant"
      ) %>%
      # First calculate copies per replicate using standard curve
      mutate(copy_est = 10^((Ct_Value - b) / m)) %>%
      group_by(Sample_Name) %>%
      summarise(
        # Calculate mean copies across replicates (not sum)
        mean_copies = mean(copy_est, na.rm = TRUE),
        sd_copies = sd(copy_est, na.rm = TRUE),
        # Keep Ct values for reference but don't use for abundance
        mean_Ct = mean(Ct_Value, na.rm = TRUE), 
        sd_Ct = sd(Ct_Value, na.rm = TRUE),
        n_reps = n(),
        .groups = "drop"
      )
      cleanup_tempfiles(temp_dir)
  } else {
    scale = NA
  }

  scale = scale %>% mutate(log2_qPCR_mean = ifelse(mean_copies > 0, log2(mean_copies), NA)) %>%
                    mutate(log2_qPCR_sd = ifelse(mean_copies > 0, log2(sd_copies), NA)) %>% 
                    mutate(log10_qPCR_mean = ifelse(mean_copies > 0, log10(mean_copies), NA)) %>%
                    mutate(log10_qPCR_sd = ifelse(mean_copies > 0, log10(sd_copies), NA))
  # Metadata
  meta_csv     <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_zip, meta_csv)
  metadata  <- read.csv(meta_con) %>% as.data.frame() %>% rename(Accession = Run)

  metadata <- metadata %>%
    mutate(Sample_name = str_extract(Library.Name, "EH\\d{2}"))

  scale <- scale %>%
    mutate(Sample_name = str_replace_all(str_extract(Sample_Name, "EH \\d{2}"), " ", ""))
  
  scale <- scale %>%
  left_join(
    metadata %>% select(Sample_name, Accession),
    by = "Sample_name"
  )

  # ----- Original counts from CSV.zip -----
  if (file.exists(orig_counts_zip)) {
        orig_csv <- unzip(orig_counts_zip, list = TRUE)$Name[1]
        orig_con <- unz(orig_counts_zip, orig_csv)
        orig_mat <- read.csv(orig_con, row.names = 1, check.names = FALSE)
        counts_original <- as.data.frame(orig_mat)
        counts_original$ASV <- rownames(counts_original)
        rownames(counts_original) = counts_original$Sequence
        
        # Taxa
        taxonomy_cols <- c("ASV", "Sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Haplotype")
        tax_original = counts_original[, names(counts_original) %in% taxonomy_cols]
        tax_original$species = tax_original$species
        tax_original$Species = NULL
        tax_original = make_taxa_label(tax_original)
        tax_original$Species = tax_original$species
        tax_original$species = NULL
        
        # Counts
        drop <- c("ASV", "Sequence", "Kingdom", "Phylum", "Class", 
                "Order", "Family", "Genus", "Species", "Haplotype")
        counts_original <- counts_original %>% select(-all_of(drop))
        counts_original <- as.data.frame(t(counts_original))

        counts_original <- counts_original %>%
            tibble::rownames_to_column("rowname") %>%
            mutate(
                Sample_name = case_when(
                    str_detect(rowname, "EH\\d{2}") ~ str_extract(rowname, "EH\\d{2}"),
                    str_detect(rowname, "MD1b") ~ "MD1b",
                    TRUE ~ NA_character_)
            ) %>%
            select(-rowname) %>%
            column_to_rownames("Sample_name")

                # ----- rdp16 -----
        if (!file.exists(file.path(local,"rdp16classified_original.csv.zip"))) {
          if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
            required_pkgs <- c("dada2", "Biostrings")
            missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
            if (length(missing_pkgs) > 0) {
              stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                  ". Please install them before running this function.")
            }
            seqs <- Biostrings::DNAStringSet(colnames(counts_original))
            rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
            tax_original2 = make_taxa_label(rdpclassified) 
            write.csv(tax_original2, file = file.path(local, "rdp16classified_original.csv"), row.names = TRUE)
          } else {
            stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
          }
        } else {
          tax_original2 <- read_zipped_table(file.path(local, "rdp16classified_original.csv.zip"), row.names = 1)
        }

        if (!raw) {
            aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
            counts_original = aligned$counts_original
            counts_original2 = aligned$counts_original
            matched_taxa <- tax_original$Taxa[match(colnames(counts_original), rownames(tax_original))]
            matched_taxa2 <- tax_original2$Taxa[match(colnames(counts_original2), rownames(tax_original2))]
            colnames(counts_original) <- matched_taxa
            colnames(counts_original2) <- matched_taxa2
            counts_original <- collapse_duplicate_columns_exact(counts_original)
            counts_original2 <- collapse_duplicate_columns_exact(counts_original2)
            original_names <- colnames(counts_original)
            original_names2 <- colnames(counts_original2)
            counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
            counts_original2 <- as.data.frame(lapply(counts_original2, as.numeric), row.names = rownames(counts_original2), col.names = original_names2, check.names = FALSE)
            proportions_original2 <- sweep(counts_original2, 1, rowSums(counts_original2), '/')
        }

        # Calculate proportions
        proportions_original <- sweep(counts_original, MARGIN = 1, STATS = rowSums(counts_original), FUN = "/")
    } else {
        counts_original = read_zipped_table(deletewhenfinished) %>% as.data.frame() 
        if (!raw) {
          aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
          counts_original = aligned$counts_original
          original_names <- colnames(counts_original)
          counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
        }
        proportions_original = sweep(counts_original, 1, rowSums(counts_original), "/")
    }

  # ----- Reprocessed counts from RDS ZIP -----
  if (file.exists(repro_counts_rds_zip)) {
    temp_dir <- tempfile("repro")
    dir.create(temp_dir)
    unzipped <- unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
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
    unzipped <- unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))
    tax_reprocessed = make_taxa_label(tax_reprocessed)
    
    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
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
    cleanup_tempfiles(temp_dir)

    if (!raw) {
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
        proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    }
  } else {
    counts_reprocessed = NA
    proportions_reprocessed = NA
    tax_reprocessed = NA
    counts_reprocessed2 = NA
    proportions_reprocessed2 = NA
    tax_reprocessed2 = NA
  }

  if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      proportions_original = fill_na_zero_numeric(proportions_original)
      counts_original2 = fill_na_zero_numeric(counts_original2)
      proportions_original2 = fill_na_zero_numeric(proportions_original2)
  }

  # ----- Return structured list -----
  return(list(
      counts = list(
          original = list(rdp19 = counts_original, rdp16 = counts_original2),
          reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
      ),
      proportions = list(
          original = list(rdp19 = proportions_original, rdp16 = proportions_original2),
          reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
      ),
      tax = list(
          original = list(rdp19 = tax_original, rdp16 = tax_original2),
          reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
      ),
      scale = scale,
      metadata = metadata
  ))
}
