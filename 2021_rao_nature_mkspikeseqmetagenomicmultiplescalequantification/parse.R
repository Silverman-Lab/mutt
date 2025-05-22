<<<<<<< Updated upstream
parse_2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification <- function(raw = FALSE, originaltax = 'qiime', align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
=======
parse_2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification <- function() {
  
  # -------------------------------------------------------------
  # 1) Check and load required packages
  # -------------------------------------------------------------
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
>>>>>>> Stashed changes
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
<<<<<<< Updated upstream
  if (!is.logical(raw)) {
    stop("raw must be a logical value")
  }
  if (!is.logical(align)) {
    stop("align must be a logical value")
  }
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification")

  # ----- File paths -----
  counts_zip          <- file.path(local, "combined_p1_p2.xlsx.zip")
  supplemental_zip    <- file.path(local, "41586_2021_3241_MOESM4_ESM.xlsx.zip")
  sex_delivery_zip    <- file.path(local, "SI_data3_sex_and_delivery_data.csv.zip")
  diet_data_zip       <- file.path(local, "SI_data2_diet_data.csv.zip")
  meds_data_zip       <- file.path(local, "SI_data1_allMeds_jan2020.xlsx.zip")
  sra_metadata_zip    <- file.path(local, "SraRunTable (38).csv.zip")
  repro_counts_rds_zip<- file.path(local, "PRJEB36435_dada2_counts.rds.zip")
  repro_tax_zip       <- file.path(local, "PRJEB36435_dada2_taxa.rds.zip")

  safe_read_zip <- function(zip_path, is_xlsx = FALSE, sheet = 1) {
    if (!file.exists(zip_path)) return(NULL)
    zfiles <- unzip(zip_path, list = TRUE)$Name
    if (length(zfiles) < 1) return(NULL)
    file_to_read <- zfiles[1]
    td <- tempdir()
    unzip(zip_path, files = file_to_read, exdir = td, overwrite = TRUE)
    fpath <- file.path(td, file_to_read)
    if (!file.exists(fpath)) return(NULL)
    df <- tryCatch(
      {
        if (is_xlsx) {
          read_excel(fpath, sheet = sheet) %>%
            as.data.frame()
        } else {
          read.csv(fpath, stringsAsFactors = FALSE) %>%
            as.data.frame()
        }
      },
      error = function(e) NULL
    )
    return(df)
  }
  make_proportions <- function(count_df) {
    if (is.null(count_df) || nrow(count_df) == 0) return(NULL)

    row_sums <- rowSums(count_df, na.rm = TRUE)
    prop_df <- sweep(count_df, 1, row_sums, "/")
    return(prop_df)
  }
    # ---- Taxonomy table constructor ----
  make_mock_tax_table <- function(df) {
    df <- df %>%
      filter(!taxonomy %in% c("Observed_total", "Expected_total")) %>%
      mutate(
        Genus   = sub("\\..*", "", taxonomy),
        Species = sub(".*?\\.\\s*", "", taxonomy),
        Kingdom = NA, Phylum = NA, Class = NA, Order = NA, Family = NA
      ) %>%
      mutate(across(everything(), ~ ifelse(. == "", NA, .))) %>%
      mutate(Taxa = paste0("g_", Genus, "_", Species)) %>%
      select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Taxa)
=======
  
  library(tidyverse)
  library(readxl)
  
  # -------------------------------------------------------------
  # 2) Read Counts + Taxonomy from combined_p1_p2.xlsx.zip
  # -------------------------------------------------------------
  counts_zip <- "2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification/combined_p1_p2.xlsx.zip"
  
  counts_bacteria <- NULL
  counts_fungi    <- NULL
  taxa_bacteria   <- NULL
  taxa_fungi      <- NULL
  
  if (!file.exists(counts_zip)) {
    warning("Counts file not found: ", counts_zip)
  } else {
    tmp_dir <- tempdir()
    unzip(counts_zip, exdir = tmp_dir)
    excel_file <- list.files(tmp_dir, pattern = "\\.xlsx$", full.names = TRUE)
>>>>>>> Stashed changes
    
    return(df)
  }

  # ---- Process one mock OTU table (taxonomy + scale extraction + OTU formatting) ----
  process_mock_table <- function(df) {
    # Separate scale rows
    scale_df <- df %>%
      filter(taxonomy %in% c("Observed_total", "Expected_total")) %>%
      column_to_rownames("taxonomy") %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column("Sample_name") %>%
      rename(Observed_total = Observed_total, Expected_total = Expected_total) %>%
      mutate(log2_Observed_total = ifelse(Observed_total > 0, log2(Observed_total), NA),
             log2_Expected_total = ifelse(Expected_total > 0, log2(Expected_total), NA),
             log10_Observed_total = ifelse(Observed_total > 0, log10(Observed_total), NA),
             log10_Expected_total = ifelse(Expected_total > 0, log10(Expected_total), NA))

    # Process taxonomy table
    tax_df <- make_mock_tax_table(df)

    # Subset OTU table (excluding scale rows)
    otu_df <- df %>%
      filter(!taxonomy %in% c("Observed_total", "Expected_total", "log2_Observed_total", "log2_Expected_total", "log10_Observed_total", "log10_Expected_total")) %>%
      mutate(Taxa = paste0("g_", sub("\\..*", "", taxonomy), "_", sub(".*?\\.\\s*", "", taxonomy))) %>%
      select(-taxonomy) %>%
      column_to_rownames("Taxa")

    return(list(otu = otu_df, tax = tax_df, scale = scale_df))
  }

  read_otu_sheet <- function(sheetname) {
    df <- tryCatch(readxl::read_excel(file, sheet = sheetname), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    if (all(c("qiime2_sklearn_taxonomy", "confidence", "OTU_ID") %in% colnames(df))) {
      tax <- df %>%
        select(OTU_ID, qiime_sklearn = qiime2_sklearn_taxonomy, confidence)
      otu <- df %>%
        select(-qiime2_sklearn_taxonomy, -confidence)
      return(list(otu = as.data.frame(t(otu)), tax = make_taxa_label(build_taxonomy_table(tax, method = "qiime"))))
    }
    return(list(otu = df, tax = NULL))
  }

  if (!file.exists(counts_zip)) {
    warning("Counts file not found: ", counts_zip)
  } else {
    tmp_dir <- tempdir()
    unzip(counts_zip, exdir = tmp_dir)
    excel_file <- list.files(tmp_dir, pattern = "combined_p1_p2\\.xlsx$", full.names = TRUE)
    sheets <- c("fungi", "bacteria", "archaea", "archaea_side")
    data_list <- lapply(sheets, function(sh) {
      possible_data <- tryCatch(
        read_excel(excel_file, sheet = sh),
        error = function(e) NULL
      )
      return(possible_data)
    })
    names(data_list) <- sheets
    taxonomy_columns <- c(
      "qiime_sklearn", "qiime_confidence",
      "vsearch_usearchglobal", "vsearch_identity",
      "usearch_utax", "usearch_sintax", "usearch_sintax_80%"
    )
    data_fungi    <- data_list$fungi
    data_bacteria <- data_list$bacteria
    if (is.null(data_fungi))    data_fungi <- data.frame()
    if (is.null(data_bacteria)) data_bacteria <- data.frame()
    
    if (nrow(data_bacteria) > 0) {
      # Bacteria
      counts_bacteria <- data_bacteria %>%
        select(OTU_ID, everything()) %>%
<<<<<<< Updated upstream
        select(-any_of(taxonomy_columns)) 
=======
        select(-any_of(taxonomy_columns)) %>%
        column_to_rownames("OTU_ID")
      
      taxa_bacteria <- data_bacteria %>%
        select(any_of(c("OTU_ID", taxonomy_columns))) %>%
        column_to_rownames("OTU_ID")
    }
    
    if (nrow(data_fungi) > 0) {
      # Fungi
      counts_fungi <- data_fungi %>%
        select(OTU_ID, everything()) %>%
        select(-any_of(taxonomy_columns)) %>%
        column_to_rownames("OTU_ID")
      
      taxa_fungi <- data_fungi %>%
        select(any_of(c("OTU_ID", taxonomy_columns))) %>%
        column_to_rownames("OTU_ID")
    }
  }
  
  # -------------------------------------------------------------
  # 3) Read Supplemental Data from 41586_2021_3241_MOESM4_ESM.xlsx.zip
  # -------------------------------------------------------------
  supplemental_zip <- "41586_2021_3241_MOESM4_ESM.xlsx.zip"
  table_3b <- NULL
  table_3c <- NULL
  table_3d <- NULL
  scale <- NULL
  counts_bacteria_phase2 <- NULL
  counts_fungi_phase2 <- NULL
  counts_archaea <- NULL
  taxa_bacteria_phase2 <- NULL
  taxa_fungi_phase2 <- NULL
  taxa_archaea <- NULL
  
  if (file.exists(supplemental_zip)) {
    # Unzip the supplemental file
    tmp_dir <- tempdir()
    unzip(supplemental_zip, exdir = tmp_dir)
    supplemental_excel <- list.files(tmp_dir, pattern = "41586_2021_3241_MOESM4_ESM\\.xlsx$", full.names = TRUE)
    
    if (length(supplemental_excel) > 0) {
      # Sheet 4: Scale data
      scale_sheet4 <- tryCatch(
        read_excel(supplemental_excel[1], sheet = "Sheet4"),
        error = function(e) NULL
      )
      
      # Sheet 6: Scale data
      scale_sheet6 <- tryCatch(
        read_excel(supplemental_excel[1], sheet = "Sheet6"),
        error = function(e) NULL
      )
      
      # Combine scale data
      if (!is.null(scale_sheet4) && !is.null(scale_sheet6)) {
        scale <- rbind(scale_sheet4, scale_sheet6)
      } else if (!is.null(scale_sheet4)) {
        scale <- scale_sheet4
      } else if (!is.null(scale_sheet6)) {
        scale <- scale_sheet6
      }
      
      # Sheet 8: Bacteria OTU counts (1st phase)
      sheet8_data <- tryCatch(
        read_excel(supplemental_excel[1], sheet = "Sheet8"),
        error = function(e) NULL
      )
      if (!is.null(sheet8_data)) {
        counts_bacteria <- sheet8_data %>%
          select(-qiime2_sklearn_taxonomy, -qiime2_sklearn_confidence) %>%
          column_to_rownames("OTU_ID")
        
        taxa_bacteria <- sheet8_data %>%
          select(OTU_ID, qiime2_sklearn_taxonomy) %>%
          column_to_rownames("OTU_ID") %>%
          rename(taxonomy = qiime2_sklearn_taxonomy)
      }
      
      # Sheet 9: Archaea OTU counts (1st phase)
      sheet9_data <- tryCatch(
        read_excel(supplemental_excel[1], sheet = "Sheet9"),
        error = function(e) NULL
      )
      if (!is.null(sheet9_data)) {
        counts_archaea <- sheet9_data %>%
          select(-qiime2_sklearn_taxonomy, -qiime2_sklearn_confidence) %>%
          column_to_rownames("OTU_ID")
        
        taxa_archaea <- sheet9_data %>%
          select(OTU_ID, qiime2_sklearn_taxonomy) %>%
          column_to_rownames("OTU_ID") %>%
          rename(taxonomy = qiime2_sklearn_taxonomy)
      }
      
      # Sheet 10: Fungi OTU counts (1st phase)
      sheet10_data <- tryCatch(
        read_excel(supplemental_excel[1], sheet = "Sheet10"),
        error = function(e) NULL
      )
      if (!is.null(sheet10_data)) {
        counts_fungi <- sheet10_data %>%
          select(-qiime2_sklearn_taxonomy, -qiime2_sklearn_confidence) %>%
          column_to_rownames("OTU_ID")
        
        taxa_fungi <- sheet10_data %>%
          select(OTU_ID, qiime2_sklearn_taxonomy) %>%
          column_to_rownames("OTU_ID") %>%
          rename(taxonomy = qiime2_sklearn_taxonomy)
      }
      
      # Sheet 11: Bacteria OTU counts (2nd phase)
      sheet11_data <- tryCatch(
        read_excel(supplemental_excel[1], sheet = "Sheet11"),
        error = function(e) NULL
      )
      if (!is.null(sheet11_data)) {
        counts_bacteria_phase2 <- sheet11_data %>%
          select(-qiime2_sklearn_taxonomy, -qiime2_sklearn_confidence) %>%
          column_to_rownames("OTU_ID")
        
        taxa_bacteria_phase2 <- sheet11_data %>%
          select(OTU_ID, qiime2_sklearn_taxonomy) %>%
          column_to_rownames("OTU_ID") %>%
          rename(taxonomy = qiime2_sklearn_taxonomy)
      }
      
      # Sheet 12: Fungi OTU counts (2nd phase)
      sheet12_data <- tryCatch(
        read_excel(supplemental_excel[1], sheet = "Sheet12"),
        error = function(e) NULL
      )
      if (!is.null(sheet12_data)) {
        counts_fungi_phase2 <- sheet12_data %>%
          select(-qiime2_sklearn_taxonomy, -qiime2_sklearn_confidence) %>%
          column_to_rownames("OTU_ID")
        
        taxa_fungi_phase2 <- sheet12_data %>%
          select(OTU_ID, qiime2_sklearn_taxonomy) %>%
          column_to_rownames("OTU_ID") %>%
          rename(taxonomy = qiime2_sklearn_taxonomy)
      }
      
      # Subtable 3b: OTU table of bacterial mock communities
      sheet_data <- read_excel(supplemental_excel[1], sheet = "Sheet3", col_names = FALSE)
      start_row_3b <- which(sheet_data[[1]] == "3b. OTU table of bacterial mock communities") + 2
      table_3b <- read_excel(supplemental_excel[1], sheet = "Sheet3", skip = start_row_3b, n_max = 13)
      table_3b <- table_3b[1:10, ]  # First 10 rows are the taxonomy data
      colnames(table_3b)[1] <- "taxonomy"
      
      # Subtable 3c: OTU table of fungal mock communities
      start_row_3c <- which(sheet_data[[1]] == "3c. OTU table of fungal mock communities") + 2
      table_3c <- read_excel(supplemental_excel[1], sheet = "Sheet3", skip = start_row_3c, n_max = 10)
      table_3c <- table_3c[1:9, ]  # First 9 rows are the taxonomy data
      colnames(table_3c)[1] <- "taxonomy"
      
      # Subtable 3d: OTU table of archaeal mock communities
      start_row_3d <- which(sheet_data[[1]] == "3d. OTU table of archaeal mock communities") + 2
      table_3d <- read_excel(supplemental_excel[1], sheet = "Sheet3", skip = start_row_3d, n_max = 5)
      table_3d <- table_3d[1:4, ]  # First 4 rows are the taxonomy data
      colnames(table_3d)[1] <- "taxonomy"
    } else {
      warning("Supplemental Excel file not found after unzipping.")
>>>>>>> Stashed changes
    }
  } else {
    warning("Supplemental zip file not found: ", supplemental_zip)
  }
  
  # -------------------------------------------------------------
  # 4) Calculate Proportions
  # -------------------------------------------------------------
  proportions_bacteria <- NULL
  proportions_fungi <- NULL
  proportions_bacteria_phase2 <- NULL
  proportions_fungi_phase2 <- NULL
  proportions_archaea <- NULL
  
  if (!is.null(counts_bacteria) && nrow(counts_bacteria) > 0) {
    counts_bacteria_num <- counts_bacteria %>%
      mutate_all(as.numeric)
    row_sums_bac <- rowSums(counts_bacteria_num, na.rm = TRUE)
    proportions_bacteria <- sweep(counts_bacteria_num, 1, row_sums_bac, "/")
    proportions_bacteria[is.na(proportions_bacteria)] <- 0
  }
  
  if (!is.null(counts_fungi) && nrow(counts_fungi) > 0) {
    counts_fungi_num <- counts_fungi %>%
      mutate_all(as.numeric)
    row_sums_fun <- rowSums(counts_fungi_num, na.rm = TRUE)
    proportions_fungi <- sweep(counts_fungi_num, 1, row_sums_fun, "/")
    proportions_fungi[is.na(proportions_fungi)] <- 0
  }
  
  if (!is.null(counts_bacteria_phase2) && nrow(counts_bacteria_phase2) > 0) {
    counts_bacteria_phase2_num <- counts_bacteria_phase2 %>%
      mutate_all(as.numeric)
    row_sums_bac_p2 <- rowSums(counts_bacteria_phase2_num, na.rm = TRUE)
    proportions_bacteria_phase2 <- sweep(counts_bacteria_phase2_num, 1, row_sums_bac_p2, "/")
    proportions_bacteria_phase2[is.na(proportions_bacteria_phase2)] <- 0
  }
  
  if (!is.null(counts_fungi_phase2) && nrow(counts_fungi_phase2) > 0) {
    counts_fungi_phase2_num <- counts_fungi_phase2 %>%
      mutate_all(as.numeric)
    row_sums_fun_p2 <- rowSums(counts_fungi_phase2_num, na.rm = TRUE)
    proportions_fungi_phase2 <- sweep(counts_fungi_phase2_num, 1, row_sums_fun_p2, "/")
    proportions_fungi_phase2[is.na(proportions_fungi_phase2)] <- 0
  }
  
  if (!is.null(counts_archaea) && nrow(counts_archaea) > 0) {
    counts_archaea_num <- counts_archaea %>%
      mutate_all(as.numeric)
    row_sums_arc <- rowSums(counts_archaea_num, na.rm = TRUE)
    proportions_archaea <- sweep(counts_archaea_num, 1, row_sums_arc, "/")
    proportions_archaea[is.na(proportions_archaea)] <- 0
  }
  
  # -------------------------------------------------------------
  # 5) Read and merge metadata
  # -------------------------------------------------------------
  sex_delivery_zip <- "SI_data3_sex_and_delivery_data.csv.zip"
  diet_data_zip    <- "SI_data2_diet_data.csv.zip"
  meds_data_zip    <- "SI_data1_allMeds_jan2020.xlsx.zip"
  sra_metadata_zip <- "sra_metadata.csv.zip"  # New metadata file
  
  metadata <- NULL
  
  safe_read_zip <- function(zip_path, is_xlsx = FALSE, sheet = 1) {
    if (!file.exists(zip_path)) return(NULL)
    zfiles <- unzip(zip_path, list = TRUE)$Name
    if (length(zfiles) < 1) return(NULL)
    file_to_read <- zfiles[1]
    td <- tempdir()
    unzip(zip_path, files = file_to_read, exdir = td, overwrite = TRUE)
    fpath <- file.path(td, file_to_read)
    
<<<<<<< Updated upstream
    if (nrow(data_fungi) > 0) {
      # Fungi
      counts_fungi <- data_fungi %>%
        select(OTU_ID, everything()) %>%
        select(-any_of(taxonomy_columns)) 
    }
  }

  taxa_bacteria <- suppressWarnings(suppressMessages(build_taxonomy_table(data_list[["bacteria"]], method = originaltax)))
  taxa_fungi    <- suppressWarnings(suppressMessages(build_taxonomy_table(data_list[["fungi"]], method = originaltax)))

  taxa_bacteria = make_taxa_label(taxa_bacteria)
  taxa_fungi    = make_taxa_label(taxa_fungi)

  unzip(supplemental_zip, exdir = tmp_dir)
  supplemental_excel <- list.files(tmp_dir, pattern = "41586_2021_3241_MOESM4_ESM\\.xlsx$", full.names = TRUE)

  if (length(supplemental_excel) == 0) {
    warning("Supplemental Excel file not found after unzipping.")
  } else {
    file <- supplemental_excel[1]

    # Read Scale sheets
    scale_16S <- readxl::read_excel(file, sheet = "sTable4") 
    scale_ITS <- readxl::read_excel(file, sheet = "sTable6") 
    scale <- full_join(scale_16S, scale_ITS, by = "Sample_name") %>%
      mutate(across(-Sample_name, as.numeric))  %>% select(-c(`Sample description`)) %>%
      rename(
        MK_spike_univ16S = `MK-SpikeSeq-univ16S`,
        MK_spike_arch16S = `MK-SpikeSeq-arch16S`, 
        MK_spike_ITS1 = `MK-SpikeSeq-ITS1`,
        FACS_prokaryote = `FACS-prokaryote`,
        FACS_fungi_gate1 = `FACS-fungi-gate1`,
        FACS_fungi_gate2 = `FACS-fungi-gate2`, 
        FACS_fungi_gate3 = `FACS-fungi-gate3`,
        total_DNA_PicoGreen_rep1 = `total-DNA-PicoGreen-(ng/uL)-Replicate1`,
        total_DNA_PicoGreen_rep2 = `total-DNA-PicoGreen-(ng/uL)-Replicate2`,
        univ16S_Ct_mean = `rDNA-specific-qPCR-(Ct summary of 3 replicates)-univ16S-Ct-mean`,
        univ16S_Ct_SD = `rDNA-specific-qPCR-(Ct summary of 3 replicates)-univ16S-Ct-SD`,
        arch16S_Ct_mean = `rDNA-specific-qPCR-(Ct summary of 3 replicates)-arch16S-Ct-mean`,
        arch16S_Ct_SD = `rDNA-specific-qPCR-(Ct summary of 3 replicates)-arch16S-Ct-SD`,
        ITS1_Ct_mean = `rDNA-specific-qPCR-(Ct summary of 3 replicates)-ITS1-Ct-mean`,
        ITS1_Ct_SD = `rDNA-specific-qPCR-(Ct summary of 3 replicates)-ITS1-Ct-SD`,
        ITS1_qPCR_Ct = `ITS1-qPCR-Ct`,
        ITS1_total_reads = `ITS1-total-reads`,
        ITS1_spikein_reads = `ITS1-spikein-reads`,
        ITS1_total_abundance = `ITS1-total-abundance`
      ) %>%
      mutate(
        log2_MK_spike_univ16S = ifelse(MK_spike_univ16S > 0, log2(MK_spike_univ16S), NA),
        log10_MK_spike_univ16S = ifelse(MK_spike_univ16S > 0, log10(MK_spike_univ16S), NA),
        log2_MK_spike_arch16S = ifelse(MK_spike_arch16S > 0, log2(MK_spike_arch16S), NA),
        log10_MK_spike_arch16S = ifelse(MK_spike_arch16S > 0, log10(MK_spike_arch16S), NA),
        log2_MK_spike_ITS1 = ifelse(MK_spike_ITS1 > 0, log2(MK_spike_ITS1), NA),
        log10_MK_spike_ITS1 = ifelse(MK_spike_ITS1 > 0, log10(MK_spike_ITS1), NA),
        log2_FACS_prokaryote = ifelse(FACS_prokaryote > 0, log2(FACS_prokaryote), NA),
        log10_FACS_prokaryote = ifelse(FACS_prokaryote > 0, log10(FACS_prokaryote), NA),
        log2_FACS_fungi_gate1 = ifelse(FACS_fungi_gate1 > 0, log2(FACS_fungi_gate1), NA),
        log10_FACS_fungi_gate1 = ifelse(FACS_fungi_gate1 > 0, log10(FACS_fungi_gate1), NA),
        log2_FACS_fungi_gate2 = ifelse(FACS_fungi_gate2 > 0, log2(FACS_fungi_gate2), NA),
        log10_FACS_fungi_gate2 = ifelse(FACS_fungi_gate2 > 0, log10(FACS_fungi_gate2), NA),
        log2_FACS_fungi_gate3 = ifelse(FACS_fungi_gate3 > 0, log2(FACS_fungi_gate3), NA),
        log10_FACS_fungi_gate3 = ifelse(FACS_fungi_gate3 > 0, log10(FACS_fungi_gate3), NA),
        log2_total_DNA_PicoGreen_rep1 = ifelse(total_DNA_PicoGreen_rep1 > 0, log2(total_DNA_PicoGreen_rep1), NA),
        log10_total_DNA_PicoGreen_rep1 = ifelse(total_DNA_PicoGreen_rep1 > 0, log10(total_DNA_PicoGreen_rep1), NA),
        log2_total_DNA_PicoGreen_rep2 = ifelse(total_DNA_PicoGreen_rep2 > 0, log2(total_DNA_PicoGreen_rep2), NA),
        log10_total_DNA_PicoGreen_rep2 = ifelse(total_DNA_PicoGreen_rep2 > 0, log10(total_DNA_PicoGreen_rep2), NA),
        log2_ITS1_total_abundance = ifelse(ITS1_total_abundance > 0, log2(ITS1_total_abundance), NA),
        log10_ITS1_total_abundance = ifelse(ITS1_total_abundance > 0, log10(ITS1_total_abundance), NA)
      )

   # Read OTU sheets
    otu_data <- list(
      originalcounts_bacteria        =read_otu_sheet("sTable8"),
      originalcounts_archaea         =read_otu_sheet("sTable9"),
      originalcounts_fungi           =read_otu_sheet("sTable10"),
      originalcounts_bacteria_phase2 =read_otu_sheet("sTable11"),
      originalcounts_fungi_phase2    =read_otu_sheet("sTable12")
    )

    otu_tables     <- map(otu_data, "otu")
    otu_tables     <- map(otu_tables, first_row_to_colnames)
    taxonomy_tables <- map(otu_data, "tax")

    # ---- Load Sheet3 and extract mock_metadata ----
    sheet3 <- read_excel(file, sheet = "sTable3", col_names = FALSE)
    mock_metadata <- sheet3 %>%
      slice(1:40) %>%
      rename_with(~ make.names(., unique = TRUE)) %>%
      mutate(across(everything(), readr::parse_guess))
    colnames(mock_metadata) <- as.character(unlist(mock_metadata[1, ])) 
    mock_metadata <- mock_metadata[-1, ] %>% select(1:4)
    mock_metadata <- mock_metadata %>% mutate(across(everything(), ~ na_if(., "/")))

    # ---- Function to read each OTU table from its labeled section ----
    read_mock_table <- function(label, n_max) {
      start <- which(sheet3[[1]] == label) 
      tbl <- read_excel(file, sheet = "sTable3", skip = start, n_max = n_max)
      tbl <- tbl[1:(n_max - 1), ]
      return(tbl)
=======
    if (!file.exists(fpath)) return(NULL)
    
    df <- tryCatch(
      {
        if (is_xlsx) {
          read_excel(fpath, sheet = sheet) %>%
            as.data.frame()
        } else {
          read.csv(fpath, stringsAsFactors = FALSE) %>%
            as.data.frame()
        }
      },
      error = function(e) NULL
    )
    return(df)
  }
  
  # Read existing metadata files
  sex_delivery_data <- safe_read_zip(sex_delivery_zip, is_xlsx = FALSE)
  diet_data         <- safe_read_zip(diet_data_zip, is_xlsx = FALSE)
  meds_data         <- safe_read_zip(meds_data_zip, is_xlsx = TRUE, sheet = 1)
  
  # Read new SRA metadata
  sra_metadata <- safe_read_zip(sra_metadata_zip, is_xlsx = FALSE)
  
  # Rename baby_id to id in sex_delivery_data if it exists
  if (!is.null(sex_delivery_data)) {
    if ("baby_id" %in% colnames(sex_delivery_data)) {
      sex_delivery_data <- sex_delivery_data %>%
        rename(id = baby_id)
>>>>>>> Stashed changes
    }
  }
  
  # Merge existing metadata
  if (!is.null(sex_delivery_data)) {
    metadata <- sex_delivery_data
    if (!is.null(diet_data)) {
      metadata <- metadata %>%
        full_join(diet_data, by = "id")
    }
    if (!is.null(meds_data)) {
      metadata <- metadata %>%
        full_join(meds_data, by = "id")
    }
  }
  
  # Note: The sample_name column in sra_metadata (e.g., "Ca-10-1-R1-low_S70_L001_R1") may not directly match the 'id' column in the existing metadata.
  # may need to clean and standardize the sample_name and id columns to ensure proper merging (e.g., by extracting a common identifier or removing suffixes).
  
  # Merge SRA metadata using sample_name
  if (!is.null(sra_metadata) && !is.null(metadata)) {
    # Attempt to merge using sample_name and id
    metadata <- metadata %>%
      full_join(sra_metadata, by = c("id" = "sample_name"))
  } else if (!is.null(sra_metadata) && is.null(metadata)) {
    # If no existing metadata, use SRA metadata as the base
    metadata <- sra_metadata
  }
  
  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata)
  }

<<<<<<< Updated upstream
    # ---- Load the bacterial, fungal, and fecal test OTU tables ----
    mock_otu_tables <- list(
      originalcounts_mockbacteria = read_mock_table("3b. OTU table of bacterial mock communities", 13),
      originalcounts_mockfungi    = read_mock_table("3c. OTU table of fungal mock communities", 13),
      originalcounts_mockfecal  = read_mock_table("3d. OTU table of test fecal samples", 107)
    )

    mock_otu_tables$originalcounts_mockfecal <- mock_otu_tables$originalcounts_mockfecal %>%
      separate("taxonomy-assignment",
              into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
              sep = ";", fill = "right") %>%
      mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .))) 

    originaltax_mockfecal <- mock_otu_tables$originalcounts_mockfecal %>%
      select(Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
      as.data.frame() %>%
      make_taxa_label()

    taxa <- originaltax_mockfecal$Taxa
    mock_otu_tables$originalcounts_mockfecal$Taxa <- taxa
    mock_otu_tables$originalcounts_mockfecal <- mock_otu_tables$originalcounts_mockfecal %>%
      group_by(Taxa) %>%
      summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
      as.data.frame()
    mock_otu_tables$originalcounts_mockfecal <- mock_otu_tables$originalcounts_mockfecal %>%
      column_to_rownames("Taxa")

    # ---- Apply to both mockbacteria and mockfungi ----
    mock_bacteria_result <- process_mock_table(mock_otu_tables$originalcounts_mockbacteria)
    mock_fungi_result    <- process_mock_table(mock_otu_tables$originalcounts_mockfungi)
    mock_otu_tables$originalcounts_mockbacteria <- as.data.frame(t(mock_bacteria_result$otu))
    mock_otu_tables$originalcounts_mockfungi    <- as.data.frame(t(mock_fungi_result$otu))
    originaltax_mockbacteria <- mock_bacteria_result$tax
    originaltax_mockfungi    <- mock_fungi_result$tax
    scale_mockbacteria <- mock_bacteria_result$scale
    scale_mockfungi    <- mock_fungi_result$scale
    scale_mockfecal    <- NA
  }
  
  # Read existing metadata files
  sex_delivery_data <- safe_read_zip(sex_delivery_zip, is_xlsx = FALSE)
  diet_data         <- safe_read_zip(diet_data_zip, is_xlsx = FALSE)
  meds_data         <- safe_read_zip(meds_data_zip, is_xlsx = TRUE, sheet = 1)
  sra_metadata      <- safe_read_zip(sra_metadata_zip, is_xlsx = FALSE)
  

  # I dont know if we can merge all of this?
  # Note: The sample_name column in sra_metadata (e.g., "Ca-10-1-R1-low_S70_L001_R1") may not directly match the 'id' column in the existing metadata.
  # may need to clean and standardize the sample_name and id columns to ensure proper merging (e.g., by extracting a common identifier or removing suffixes).
  # I Dont think we can link the accessions to the metadata.

  # if (!is.null(sex_delivery_data)) {
  #   if ("baby_id" %in% colnames(sex_delivery_data)) {
  #     sex_delivery_data <- sex_delivery_data %>%
  #       rename(id = baby_id)
  #   }
  #   metadata <- sex_delivery_data
  #   if (!is.null(diet_data)) {
  #     metadata <- metadata %>%
  #       full_join(diet_data, by = "id")
  #   }
  #   if (!is.null(meds_data)) {
  #     metadata <- metadata %>%
  #       full_join(meds_data, by = "id")
  #   }
  # }
  
  # # Merge SRA metadata using sample_name
  # if (!is.null(sra_metadata) && !is.null(metadata)) {
  #   # Attempt to merge using sample_name and id
  #   metadata <- metadata %>%
  #     full_join(sra_metadata, by = c("id" = "sample_name"))
  #   metadata <- as.data.frame(metadata)
  # }

  metadata = sra_metadata %>%
      as.data.frame() %>%
      rename(Accession = Run) %>%
  mutate(
    sequencer_type = case_when(
      str_detect(Instrument, regex("miseq", ignore_case = TRUE)) ~ "Miseq",
      str_detect(Instrument, regex("nextseq", ignore_case = TRUE)) ~ "Nextseq",
      TRUE ~ "Other"
    )
  ) %>%
  mutate(
    sample_type = case_when(
      str_detect(Sample_name, "_bac16SV3V4") ~ "bac16SV3V4",
      str_detect(Sample_name, "_ITS1") ~ "ITS1",
      str_detect(Sample_name, "_arch16S") ~ "arch16S",
      TRUE ~ NA_character_
    ),
    sample_prefix = str_remove(Sample_name, "_bac16SV3V4.*|_ITS1.*|_arch16S.*"),
    combinedphasesamplename = case_when(
      sequencer_type %in% c("Miseq", "Nextseq") ~ paste(sequencer_type, sample_prefix, sep = "_"),
      TRUE ~ NA_character_
    )
  )

  # Merge SRA metadata with scale data using sample_name
  if (!is.null(metadata) && !is.null(scale)) {
    scale <- scale %>%
      left_join(metadata %>% select(Sample_name, Accession), by = "Sample_name") %>%
      as.data.frame()
  }

  scale_expanded <- metadata %>%
  select(Accession, Sample_name, combinedphasesamplename) %>%
  mutate(
    Sample_scale = map_chr(
      Sample_name,
      ~ {
        hits <- scale$Sample_name[
          str_detect(.x, fixed(scale$Sample_name))
        ]
        if (length(hits)>=1) hits[1] else NA_character_
      }
    )
  ) %>%
  left_join(
    scale,
    by = c("Sample_scale" = "Sample_name")
  )

  scale_matched_only <- scale_expanded %>%
    filter(!is.na(Sample_scale)) 


  counts_bacteria <- counts_bacteria %>%
    left_join(taxa_bacteria %>% select(OTU_ID, Taxa), by = "OTU_ID") %>%
    select(-OTU_ID) %>%
    group_by(Taxa) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("Taxa") %>% t() %>% as.data.frame()

  counts_fungi <- counts_fungi %>%
    left_join(taxa_fungi %>% select(OTU_ID, Taxa), by = "OTU_ID") %>%
    select(-OTU_ID) %>%
    group_by(Taxa) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("Taxa") %>% t() %>% as.data.frame()

  otu_tables$originalcounts_bacteria = otu_tables$originalcounts_bacteria %>% t() %>% as.data.frame() %>% rownames_to_column("OTU_ID")
  otu_tables$originalcounts_bacteria <- otu_tables$originalcounts_bacteria %>%
    left_join(taxonomy_tables$originalcounts_bacteria %>% select(OTU_ID, Taxa), by = "OTU_ID") %>%
    select(-OTU_ID) %>%
    group_by(Taxa) %>%
    mutate(across(everything(), as.numeric)) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("Taxa") %>% t() %>% as.data.frame()

  otu_tables$originalcounts_fungi = otu_tables$originalcounts_fungi %>% t() %>% as.data.frame() %>% rownames_to_column("OTU_ID")
  otu_tables$originalcounts_fungi <- otu_tables$originalcounts_fungi %>%
    left_join(taxonomy_tables$originalcounts_fungi %>% select(OTU_ID, Taxa), by = "OTU_ID") %>%
    select(-OTU_ID) %>%
    group_by(Taxa) %>%
    mutate(across(everything(), as.numeric)) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("Taxa") %>% t() %>% as.data.frame()
    
  otu_tables$originalcounts_archaea = otu_tables$originalcounts_archaea %>% t() %>% as.data.frame() %>% rownames_to_column("OTU_ID")
  otu_tables$originalcounts_archaea <- otu_tables$originalcounts_archaea %>%
    left_join(taxonomy_tables$originalcounts_archaea %>% select(OTU_ID, Taxa), by = "OTU_ID") %>%
    select(-OTU_ID) %>%
    group_by(Taxa) %>%
    mutate(across(everything(), as.numeric)) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("Taxa") %>% t() %>% as.data.frame()
    
  otu_tables$originalcounts_bacteria_phase2 = otu_tables$originalcounts_bacteria_phase2 %>% t() %>% as.data.frame() %>% rownames_to_column("OTU_ID")
  otu_tables$originalcounts_bacteria_phase2 <- otu_tables$originalcounts_bacteria_phase2 %>%
    left_join(taxonomy_tables$originalcounts_bacteria_phase2 %>% select(OTU_ID, Taxa), by = "OTU_ID") %>%
    select(-OTU_ID) %>%
    group_by(Taxa) %>%
    mutate(across(everything(), as.numeric)) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("Taxa") %>% t() %>% as.data.frame()
  
  otu_tables$originalcounts_fungi_phase2 = otu_tables$originalcounts_fungi_phase2 %>% t() %>% as.data.frame() %>% rownames_to_column("OTU_ID")
  otu_tables$originalcounts_fungi_phase2 <- otu_tables$originalcounts_fungi_phase2 %>%
    left_join(taxonomy_tables$originalcounts_fungi_phase2 %>% select(OTU_ID, Taxa), by = "OTU_ID") %>%
    select(-OTU_ID) %>%
    group_by(Taxa) %>%
    mutate(across(everything(), as.numeric)) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("Taxa") %>% t() %>% as.data.frame()



    if (!raw) {
      # Align original bacteria counts
      aligned_bacteria = rename_and_align(counts_original = otu_tables$originalcounts_bacteria, 
                                        metadata = metadata, 
                                        scale = scale, 
                                        by_col = "Sample_name",
                                        align = align,
                                        study_name = basename(local))
      original_names <- colnames(aligned_bacteria$counts_original)
      aligned_bacteria$counts_original <- as.data.frame(lapply(aligned_bacteria$counts_original, as.numeric), row.names = rownames(aligned_bacteria$counts_original), col.names = original_names, check.names = FALSE)
      otu_tables$originalcounts_bacteria = fill_na_zero_numeric(aligned_bacteria$counts_original)
      
      # Align original archaea counts
      aligned_archaea = rename_and_align(counts_original = otu_tables$originalcounts_archaea,
                                       metadata = metadata,
                                       scale = scale,
                                       by_col = "Sample_name", 
                                       align = align,
                                       study_name = basename(local))
      original_names <- colnames(aligned_archaea$counts_original)
      aligned_archaea$counts_original <- as.data.frame(lapply(aligned_archaea$counts_original, as.numeric), row.names = rownames(aligned_archaea$counts_original), col.names = original_names, check.names = FALSE)
      otu_tables$originalcounts_archaea = fill_na_zero_numeric(aligned_archaea$counts_original)
      
      # Align original fungi counts
      aligned_fungi = rename_and_align(counts_original = otu_tables$originalcounts_fungi,
                                     metadata = metadata,
                                     scale = scale,
                                     by_col = "Sample_name",
                                     align = align,
                                     study_name = basename(local))
      original_names <- colnames(aligned_fungi$counts_original)
      aligned_fungi$counts_original <- as.data.frame(lapply(aligned_fungi$counts_original, as.numeric), row.names = rownames(aligned_fungi$counts_original), col.names = original_names, check.names = FALSE)
      otu_tables$originalcounts_fungi = fill_na_zero_numeric(aligned_fungi$counts_original)
      
      # Align phase 2 bacteria counts
      aligned_bacteria_p2 = rename_and_align(counts_original = otu_tables$originalcounts_bacteria_phase2,
                                           metadata = metadata,
                                           scale = scale,
                                           by_col = "Sample_name",
                                           align = align,
                                           study_name = basename(local))
      original_names <- colnames(aligned_bacteria_p2$counts_original)
      aligned_bacteria_p2$counts_original <- as.data.frame(lapply(aligned_bacteria_p2$counts_original, as.numeric), row.names = rownames(aligned_bacteria_p2$counts_original), col.names = original_names, check.names = FALSE)
      otu_tables$originalcounts_bacteria_phase2 = fill_na_zero_numeric(aligned_bacteria_p2$counts_original)
      
      # Align phase 2 fungi counts
      aligned_fungi_p2 = rename_and_align(counts_original = otu_tables$originalcounts_fungi_phase2,
                                        metadata = metadata,
                                        scale = scale,
                                        by_col = "Sample_name",
                                        align = align,
                                        study_name = basename(local))
      original_names <- colnames(aligned_fungi_p2$counts_original)
      aligned_fungi_p2$counts_original <- as.data.frame(lapply(aligned_fungi_p2$counts_original, as.numeric), row.names = rownames(aligned_fungi_p2$counts_original), col.names = original_names, check.names = FALSE)
      otu_tables$originalcounts_fungi_phase2 = fill_na_zero_numeric(aligned_fungi_p2$counts_original)
      
      # Align mock bacteria counts
      aligned_mock_bacteria = rename_and_align(counts_original = mock_otu_tables$originalcounts_mockbacteria,
                                             metadata = mock_metadata,
                                             scale = scale_mockbacteria,
                                             by_col = "Sample_name",
                                             align = align,
                                             study_name = basename(local))
      original_names <- colnames(aligned_mock_bacteria$counts_original)
      aligned_mock_bacteria$counts_original <- as.data.frame(lapply(aligned_mock_bacteria$counts_original, as.numeric), row.names = rownames(aligned_mock_bacteria$counts_original), col.names = original_names, check.names = FALSE)
      mock_otu_tables$originalcounts_mockbacteria = fill_na_zero_numeric(aligned_mock_bacteria$counts_original)
      
      # Align mock fungi counts
      aligned_mock_fungi = rename_and_align(counts_original = mock_otu_tables$originalcounts_mockfungi,
                                          metadata = mock_metadata,
                                          scale = scale_mockfungi,
                                          by_col = "Sample_name",
                                          align = align,
                                          study_name = basename(local))
      original_names <- colnames(aligned_mock_fungi$counts_original)
      aligned_mock_fungi$counts_original <- as.data.frame(lapply(aligned_mock_fungi$counts_original, as.numeric), row.names = rownames(aligned_mock_fungi$counts_original), col.names = original_names, check.names = FALSE)
      mock_otu_tables$originalcounts_mockfungi = fill_na_zero_numeric(aligned_mock_fungi$counts_original)
      
      # Align mock fecal counts
      aligned_mock_fecal = rename_and_align(counts_original = mock_otu_tables$originalcounts_mockfecal,
                                          metadata = mock_metadata,
                                          scale = scale_mockfecal,
                                          by_col = "Sample_name",
                                          align = align,
                                          study_name = basename(local))
      original_names <- colnames(aligned_mock_fecal$counts_original)
      aligned_mock_fecal$counts_original <- as.data.frame(lapply(aligned_mock_fecal$counts_original, as.numeric), row.names = rownames(aligned_mock_fecal$counts_original), col.names = original_names, check.names = FALSE)
      mock_otu_tables$originalcounts_mockfecal = fill_na_zero_numeric(aligned_mock_fecal$counts_original)
      
      # Align reprocessed bacteria and fungi counts
      aligned_bacteria_repro = rename_and_align(counts_reprocessed = counts_bacteria,
                                              metadata = metadata,
                                              scale = scale_matched_only,
                                              by_col = "combinedphasesamplename",
                                              align = align,
                                              study_name = basename(local))
      original_names <- colnames(aligned_bacteria_repro$reprocessed)
      aligned_bacteria_repro$reprocessed <- as.data.frame(lapply(aligned_bacteria_repro$reprocessed, as.numeric), row.names = rownames(aligned_bacteria_repro$reprocessed), col.names = original_names, check.names = FALSE)
      counts_bacteria = fill_na_zero_numeric(aligned_bacteria_repro$reprocessed)
      
      aligned_fungi_repro = rename_and_align(counts_reprocessed = counts_fungi,
                                           metadata = metadata,
                                           scale = scale_matched_only,
                                           by_col = "combinedphasesamplename",
                                           align = align,
                                           study_name = basename(local))
      original_names <- colnames(aligned_fungi_repro$reprocessed)
      aligned_fungi_repro$reprocessed <- as.data.frame(lapply(aligned_fungi_repro$reprocessed, as.numeric), row.names = rownames(aligned_fungi_repro$reprocessed), col.names = original_names, check.names = FALSE)
      counts_fungi = fill_na_zero_numeric(aligned_fungi_repro$reprocessed)
  }

  proportions_bacteria                = make_proportions(counts_bacteria)
  proportions_fungi                   = make_proportions(counts_fungi)
  originalproportions_bacteria        = make_proportions(otu_tables$originalcounts_bacteria)
  originalproportions_archaea         = make_proportions(otu_tables$originalcounts_archaea)
  originalproportions_fungi           = make_proportions(otu_tables$originalcounts_fungi)
  originalproportions_bacteria_phase2 = make_proportions(otu_tables$originalcounts_bacteria_phase2)
  originalproportions_fungi_phase2    = make_proportions(otu_tables$originalcounts_fungi_phase2)
  originalproportions_mockbacteria    = make_proportions(mock_otu_tables$originalcounts_mockbacteria)
  originalproportions_mockfungi       = make_proportions(mock_otu_tables$originalcounts_mockfungi)
  originalproportions_mockfecal       = make_proportions(mock_otu_tables$originalcounts_mockfecal)

  counts_reprocessed2 = NA
  proportions_reprocessed2 = NA
  tax_reprocessed2 = NA

  # NEED TO SEPARATE THE DIFFERENT SAMPLES INTO SPECIFIC DATASETS

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
    temp_dir <- tempfile("repro")
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
    
    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale_matched_only, by_col="Sample_name", align = align, study_name=basename(local))
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
    proportions_reprocessed = sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), FUN = "/")
    cleanup_tempfiles(temp_dir)
  }

  if (!raw) {
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
      proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
  }
  
  # -------------------------------------------------------------
  return(list(
    counts = list(
      original = list(
        combinedphases = list(
          bacteria = counts_bacteria,
          archaea = NA,
          fungi = counts_fungi
        ),
        phase1 = list(
          bacteria = otu_tables$originalcounts_bacteria,
          archaea = otu_tables$originalcounts_archaea,
          fungi = otu_tables$originalcounts_fungi
        ),
        phase2 = list(
          bacteria = otu_tables$originalcounts_bacteria_phase2,
          archaea = NA,
          fungi = otu_tables$originalcounts_fungi_phase2
        ),
        mock = list(
          bacteria = mock_otu_tables$originalcounts_mockbacteria,
          archaea = NA,
          fungi = mock_otu_tables$originalcounts_mockfungi,
          fecal = mock_otu_tables$originalcounts_mockfecal
        )
      ),
      reprocessed = list(
          rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2
      )
    ),
    proportions = list(
      original = list(
        combinedphases = list(
          bacteria = proportions_bacteria,
          archaea = NA,
          fungi = proportions_fungi
        ),
        phase1 = list(
          bacteria = originalproportions_bacteria,
          archaea = originalproportions_archaea,
          fungi = originalproportions_fungi
        ),
        phase2 = list(
          bacteria = originalproportions_bacteria_phase2,
          archaea = NA,
          fungi = originalproportions_fungi_phase2
        ),
        mock = list(
          bacteria = ifelse(is.null(mock_otu_tables$originalproportions_mockbacteria), NA, mock_otu_tables$originalproportions_mockbacteria),
          archaea = NA,
          fungi = ifelse(is.null(mock_otu_tables$originalproportions_mockfungi), NA, mock_otu_tables$originalproportions_mockfungi),
          fecal = ifelse(is.null(mock_otu_tables$originalproportions_mockfecal), NA, mock_otu_tables$originalproportions_mockfecal)
        )
      ),
      reprocessed = list(
          rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2
      )
    ),
    tax = list(
      original = list(
        combinedphases = list(
          bacteria = taxa_bacteria,
          archaea = NA,
          fungi = taxa_fungi
        ),
        phase1 = list(
          bacteria = taxonomy_tables$originaltax_bacteria,
          archaea = taxonomy_tables$originaltax_archaea,
          fungi = taxonomy_tables$originaltax_fungi
        ),
        phase2 = list(
          bacteria = taxonomy_tables$originaltax_bacteria_phase2,
          archaea = NA,
          fungi = taxonomy_tables$originaltax_fungi_phase2
        ),
        mock = list(
          bacteria = originaltax_mockbacteria,
          archaea = NA,
          fungi = originaltax_mockfungi,
          fecal = originaltax_mockfecal
        )
      ),
      reprocessed = list(
          rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2
      )
    ),
    scale = list(combinedphases = scale_matched_only,
                  mockfecal = scale_mockfecal,
                  mockfungi = scale_mockfungi,
                  mockbacteria = scale_mockbacteria),
    metadata = metadata
=======
  local               <-  "2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification/"
  repro_counts_rds_zip<-  paste0(local, "PRJEB36435_dada2_merged_nochim.rds.zip"
  repro_tax_zip       <-  paste0(local, "PRJEB36435_dada2_taxonomy_merged.rds.zip"

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
  
  # -------------------------------------------------------------
  # 6) Return final list
  # -------------------------------------------------------------
  return(list(
    counts = list(
      counts_bacteria = counts_bacteria,
      counts_fungi = counts_fungi,
      counts_archaea = counts_archaea,
      counts_bacteria_phase2 = counts_bacteria_phase2,
      counts_fungi_phase2 = counts_fungi_phase2
    ),
    proportions = list(
      proportions_bacteria = proportions_bacteria,
      proportions_fungi = proportions_fungi,
      proportions_archaea = proportions_archaea,
      proportions_bacteria_phase2 = proportions_bacteria_phase2,
      proportions_fungi_phase2 = proportions_fungi_phase2
    ),
    tax = list(
      taxa_bacteria = taxa_bacteria,
      taxa_fungi = taxa_fungi,
      taxa_archaea = taxa_archaea,
      taxa_bacteria_phase2 = taxa_bacteria_phase2,
      taxa_fungi_phase2 = taxa_fungi_phase2
    ),
    scale = scale,
    metadata = metadata,
    supplemental_tables = list(
      table_3b = table_3b,
      table_3c = table_3c,
      table_3d = table_3d
    )
>>>>>>> Stashed changes
  ))
}