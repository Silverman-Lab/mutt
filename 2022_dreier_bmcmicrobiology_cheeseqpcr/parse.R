parse_2022_dreier_bmcmicrobiology_cheeseqpcr <- function() {
  required_pkgs <- c("tibble", "tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)
  library(readxl)

    # ----- Local base directory -----
    local <- file.path("2022_dreier_bmcmicrobiology_cheeseqpcr")

    # ----- File paths -----
    metadata_zip         <- file.path(local, "SraRunTable (30).csv")
    scale_zip            <- file.path(local, "V18-22-21_htqpcr_rawdata.csv.zip")
    orig_counts_zip      <- file.path(local, "V18-22-21_ASV_counts_table.csv.zip")
    orig_prop_zip        <- NA
    repro_counts_rds_zip <- file.path(local, "PRJNA786903_dada2_merged_nochim.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA786903_dada2_taxonomy_merged.rds.zip")

  # ----- Original counts from CSV.zip -----
  if (file.exists(orig_counts_zip)) {
    orig_csv <- unzip(orig_counts_zip, list = TRUE)$Name[1]
    orig_con <- unz(orig_counts_zip, orig_csv)
    orig_mat <- read.csv(orig_con, row.names = 1, check.names = FALSE)
    counts_original <- as.data.frame(orig_mat)
    counts_original$Sequence <- rownames(counts_original)
    counts_original$ASV <- counts_original$ID
    columns_to_drop <- c("ASV", "ID", "Kingdom", "Order", "Phylum", "Family", "Genus", "Haplotype", "Class", "Species")
    # Taxa
    taxonomy_cols <- c("ASV", "Sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Haplotype")
    tax = counts_original[, names(counts_original) %in% taxonomy_cols]
    rownames(tax) <- paste0("Taxon_", seq_len(nrow(tax)))
    # Counts
    counts_original <- counts_original[, !(names(counts_original) %in% columns_to_drop)]
    counts_original = counts_original[, c("Sequence", setdiff(names(counts_original), "Sequence"))]
    rownames(counts_original) <- paste0("Taxon_", seq_len(nrow(counts_original)))
  } else {
    counts_original = NA
  }

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
  proportions_original = NA
  }

  # Scale
  if (file.exists(scale_zip)) {
    # Extract raw Fluidigm CSV from the zip archive
    scale_file <- unzip(scale_zip, list = TRUE)$Name[1]
    csv_file   <- unzip(scale_zip, files = scale_file, exdir = tempdir(), overwrite = TRUE)[1]

    # Read and flatten Fluidigm 2-row header
    raw_hdr <- read_lines(csv_file, n_max = 11)
    hdr1    <- str_split(raw_hdr[10], ",")[[1]] %>% str_trim()
    hdr2    <- str_split(raw_hdr[11], ",")[[1]] %>% str_trim()
    col_nms <- make.names(paste0(hdr1, "_", hdr2), allow_ = TRUE)

    # Load qPCR data
    df <- read_csv(csv_file, skip = 11, col_names = col_nms) 
    df <- df %>%
    rename(
        Sample_Name      = Experiment.Information_Sample...2,
        Sample_Type      = Experiment.Information_Sample...3,
        Sample_rConc_raw = Experiment.Information_Sample...4,
        Ct               = EvaGreen_Ct...7,
        Calibrated_rConc = EvaGreen_Ct...8
    ) %>%
    mutate(
        Ct               = as.numeric(Ct),
        Sample_rConc_raw = as.numeric(Sample_rConc_raw),
        Calibrated_rConc = as.numeric(Calibrated_rConc)
    )
    # 1. Total abundance from calibrated concentrations
    df %>%
        group_by(Sample_Name) %>%
        summarise(
        total_copies = sum(Calibrated_rConc, na.rm = TRUE),
        mean_copies  = mean(Calibrated_rConc, na.rm = TRUE),
        sd_copies    = sd(Calibrated_rConc, na.rm = TRUE),
        n_reps       = sum(!is.na(Calibrated_rConc)),
        .groups = "drop"
        ) 
    # 2. Fit standard curve from "Standard" wells
    stds <- df %>%
        filter(
        str_to_lower(Sample_Type) == "standard",
        !is.na(Ct),
        !is.na(Sample_rConc_raw)
        ) %>%
        mutate(log_copies = log10(Sample_rConc_raw))

    if (nrow(stds) < 2) {
        stop("✗ Not enough standard wells found to fit a Ct → copy curve.")
    }

    fit <- lm(Ct ~ log_copies, data = stds)
    m   <- coef(fit)["log_copies"]
    b   <- coef(fit)["(Intercept)"]
    message(sprintf("✓ Standard curve: Ct = %.4f × log10(copies) + %.4f", m, b))

    # 3. Apply standard curve to unknowns
    df %>%
        filter(
        !is.na(Ct),
        str_to_lower(Sample_Type) != "standard"
        ) %>%
        mutate(copy_est = 10^((Ct - b) / m)) %>%
        group_by(Sample_Name) %>%
        summarise(
        total_copies_est = sum(copy_est, na.rm = TRUE),
        mean_Ct          = mean(Ct, na.rm = TRUE),
        sd_Ct            = sd(Ct,  na.rm = TRUE),
        n_reps           = n(),
        .groups = "drop"
        ) 

    scale = df
  } else {
    scale = NA
  }

  # Metadata
  meta_csv     <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_zip, meta_csv)
  metadata  <- read.csv(meta_con) %>% as.data.frame()
  
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
  tax_reprocessed <- tax_table

  return(list(
    counts = list(
      original = counts,
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original = proportions,
      reprocessed = proportions_reprocessed
    ),
    tax = list(
      original = tax,
      reprocessed = tax_reprocessed
    ),
    scale = scale,
    metadata = metadata
  ))
}
