parse_2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification <- function() {
<<<<<<< Updated upstream
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
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
  sra_metadata_zip    <- file.path(local, "sra_metadata.csv.zip")
  repro_counts_rds_zip<- file.path(local, "PRJEB36435_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- file.path(local, "PRJEB36435_dada2_taxonomy_merged.rds.zip")

  counts_bacteria <- NULL
  counts_fungi    <- NULL
  taxa_bacteria   <- NULL
  taxa_fungi      <- NULL
  
=======
  
  # -------------------------------------------------------------
  # 1) Check and load required packages
  # -------------------------------------------------------------
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  
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
  
>>>>>>> Stashed changes
  if (!file.exists(counts_zip)) {
    warning("Counts file not found: ", counts_zip)
  } else {
    tmp_dir <- tempdir()
    unzip(counts_zip, exdir = tmp_dir)
    excel_file <- list.files(tmp_dir, pattern = "\\.xlsx$", full.names = TRUE)
    
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
  
<<<<<<< Updated upstream
=======
  # -------------------------------------------------------------
  # 3) Read Supplemental Data from 41586_2021_3241_MOESM4_ESM.xlsx.zip
  # -------------------------------------------------------------
  supplemental_zip <- "41586_2021_3241_MOESM4_ESM.xlsx.zip"
>>>>>>> Stashed changes
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
    }
  } else {
    warning("Supplemental zip file not found: ", supplemental_zip)
  }
  
<<<<<<< Updated upstream
=======
  # -------------------------------------------------------------
  # 4) Calculate Proportions
  # -------------------------------------------------------------
>>>>>>> Stashed changes
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
  
<<<<<<< Updated upstream
=======
  # -------------------------------------------------------------
  # 5) Read and merge metadata
  # -------------------------------------------------------------
  sex_delivery_zip <- "SI_data3_sex_and_delivery_data.csv.zip"
  diet_data_zip    <- "SI_data2_diet_data.csv.zip"
  meds_data_zip    <- "SI_data1_allMeds_jan2020.xlsx.zip"
  sra_metadata_zip <- "sra_metadata.csv.zip"  # New metadata file
  
>>>>>>> Stashed changes
  metadata <- NULL
  
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

>>>>>>> Stashed changes
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
  ))
}