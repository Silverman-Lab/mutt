parse_2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification <- function(raw = FALSE, originaltax = 'qiime', align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl")
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

    count_df_num <- count_df %>%
      mutate(across(everything(), as.numeric))

    row_sums <- rowSums(count_df_num, na.rm = TRUE)
    prop_df <- sweep(count_df_num, 1, row_sums, "/")
    prop_df[is.na(prop_df)] <- 0
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
        select(-any_of(taxonomy_columns)) 
    }
    
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

  tmp_dir <- tempdir()
  unzip(supplemental_zip, exdir = tmp_dir)
  supplemental_excel <- list.files(tmp_dir, pattern = "41586_2021_3241_MOESM4_ESM\\.xlsx$", full.names = TRUE)

  if (length(supplemental_excel) == 0) {
    warning("Supplemental Excel file not found after unzipping.")
  } else {
    file <- supplemental_excel[1]

    # Read Scale sheets
    scale_ <- readxl::read_excel(file, sheet = "sTable4")
    scale_ITS <- readxl::read_excel(file, sheet = "sTable6")
    scale <- bind_rows(Filter(Negate(is.null), list(scale_, scale_ITS))) %>% select(-c(`Sample description`)) %>%
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
      mutate(across(-Sample_name, as.numeric)) %>%
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
    }

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
      rename(Accession = Run)

  # Merge SRA metadata with scale data using sample_name
  if (!is.null(metadata) && !is.null(scale)) {
    scale <- scale %>%
      left_join(metadata %>% select(Sample_name, Accession), by = "Sample_name") %>%
      as.data.frame()
  }

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
                                          scale = scale_mockfungi,
                                          by_col = "Sample_name",
                                          align = align,
                                          study_name = basename(local))
      original_names <- colnames(aligned_mock_fecal$counts_original)
      aligned_mock_fecal$counts_original <- as.data.frame(lapply(aligned_mock_fecal$counts_original, as.numeric), row.names = rownames(aligned_mock_fecal$counts_original), col.names = original_names, check.names = FALSE)
      mock_otu_tables$originalcounts_mockfecal = fill_na_zero_numeric(aligned_mock_fecal$counts_original)
      
      # Align reprocessed bacteria and fungi counts
      aligned_bacteria_repro = rename_and_align(counts_reprocessed = counts_bacteria,
                                              metadata = metadata,
                                              scale = scale,
                                              by_col = "Sample_name",
                                              align = align,
                                              study_name = basename(local))
      original_names <- colnames(aligned_bacteria_repro$reprocessed)
      aligned_bacteria_repro$reprocessed <- as.data.frame(lapply(aligned_bacteria_repro$reprocessed, as.numeric), row.names = rownames(aligned_bacteria_repro$reprocessed), col.names = original_names, check.names = FALSE)
      counts_bacteria = fill_na_zero_numeric(aligned_bacteria_repro$reprocessed)
      
      aligned_fungi_repro = rename_and_align(counts_reprocessed = counts_fungi,
                                           metadata = metadata,
                                           scale = scale,
                                           by_col = "Sample_name",
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

  # NEED TO SEPARATE THE DIFFERENT SAMPLES INTO SPECIFIC DATASETS

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
    temp_dir <- tempfile("repro")
    dir.create(temp_dir)
    unzipped = unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
    counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(counts_file))

    # ----- Taxonomy reprocessed -----
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))
    
    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
      counts_reprocessed = aligned$reprocessed
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
      original_names <- colnames(counts_reprocessed)
      counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
    }

    # proportions reprocessed
    proportions_reprocessed = sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), FUN = "/")
    cleanup_tempfiles(temp_dir)
  }

  if (!raw) {
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
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
          counts = counts_reprocessed
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
          bacteria = mock_otu_tables$originalproportions_mockbacteria,
          archaea = NA,
          fungi = mock_otu_tables$originalproportions_mockfungi,
          fecal = mock_otu_tables$originalproportions_mockfecal
        )
      ),
      reprocessed = list(
          proportions = proportions_reprocessed
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
          tax = tax_reprocessed
      )
    ),
    scale = scale,
    metadata = metadata
  ))
}