parse_2020_galazzo_frontiersincellularandinfectionmicrobiology_flowandqPCRhealthy <- function() {
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  
  library(tidyverse)
  library(readxl)

  counts <- NA
  proportions <- NA
  tax <- NA
  
  # ----- Local base directory -----
  local <- file.path("2020_galazzo_frontiersincellularandinfectionmicrobiology_flowandqPCRhealthy")

  # ----- File paths -----
  metadata_zip         <- file.path(local, "SraRunTable.csv.zip")
  scale_zip            <- file.path(local, "Table 1.XLSX.zip")
  repro_counts_rds_zip <- file.path(local, "ERP108719_dada2_merged_nochim.rds.zip")
  repro_tax_zip        <- file.path(local, "ERP108719_dada2_taxonomy_merged.rds.zip")


  zip_list <- unzip(metadata_zip, list = TRUE)
  metadata_csv <- zip_list$Name[1]
  metadata_con <- unz(metadata_zip, metadata_csv)
  metadata <- read.csv(metadata_con, row.names = "person_id") %>%
    as.data.frame() %>%
    rownames_to_column("ID")
  

  zip_list <- unzip(scale_zip, list = TRUE)
  scale_xlsx <- zip_list$Name[1]  
  temp_dir <- tempdir()
  unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
  scale_path <- file.path(temp_dir, scale_xlsx)
  
  sheet3_raw <- read_xlsx(scale_path, sheet = 3)
  scale_s3_untreated <- sheet3_raw %>%
    select(
      ID = `Sample ID`,
      FACS_Mean = `Average FACS cell count/g faeces`,
      FACS_SD = `S.D.`,
      qPCR_Mean = `Copies per gram/faeces`
    ) %>%
    mutate(Treatment = "Untreated")
  
  scale_s3_pma <- sheet3_raw %>%
    select(
      ID = `Sample ID`,
      qPCR_Mean = `Copies per gram/faeces`[grep("PMA", names(sheet3_raw))],  # Adjust to match PMA column name
    ) %>%
    mutate(
      ID = paste0(ID, "-PMA"),
      FACS_Mean = NA,
      FACS_SD = NA,
      Treatment = "PMA"
    )
  
  scale_s3 <- bind_rows(scale_s3_untreated, scale_s3_pma) %>%
    select(ID, FACS_Mean, FACS_SD, qPCR_Mean) %>%
    column_to_rownames("ID")
  
  sheet8_raw <- read_xlsx(scale_path, sheet = 8)
  
  scale_s8 <- sheet8_raw %>%
    select(
      ID = `Sample ID`,
      ddPCR_Mean = `ddPCR average copies/uL DNA`,
      ddPCR_SD = `ddPCR SD copies/uL DNA`,
      qPCR_Mean = `qPCR copies/uL DNA`
    ) %>%
    column_to_rownames("ID")
  
  scale <- scale_s3 %>%
    full_join(scale_s8, by = "row.names") %>%
    column_to_rownames("Row.names")
  
  metadata <- metadata %>%
    left_join(scale_s3_untreated %>% select(ID, Treatment), by = "ID") %>%
    column_to_rownames("ID")
  
  unlink(temp_dir, recursive = TRUE)

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

  return(list(
    counts = list(
                original = counts,
                reprocessed = counts_reprocessed
    )
    proportions = list(
                original = proportions,
                reprocessed = proportions_reprocessed
    )
    tax = list(
                original = tax,
                reprocessed = tax_reprocessed
    )
    scale = scale,
    metadata = metadata
  ))
}