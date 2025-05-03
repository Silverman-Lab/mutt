parse_2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy <- function(raw = FALSE) {
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
  local <- file.path("2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy")

  # ----- File paths -----
  metadata_zip         <- file.path(local, "SraRunTable.csv.zip")
  scale_zip            <- file.path(local, "Table 1.XLSX.zip")
  repro_counts_rds_zip <- file.path(local, "ERP108719_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "ERP108719_dada2_taxa.rds.zip")


  zip_list <- unzip(metadata_zip, list = TRUE)
  metadata_csv <- zip_list$Name[1]
  metadata_con <- unz(metadata_zip, metadata_csv)
  metadata <- read.csv(metadata_con) %>%
  as.data.frame()
  if ("Sample_name" %in% names(metadata)) {
    metadata <- metadata %>%
      rename(Sample_name_existing = Sample_name)
  }
  metadata <- metadata %>%
    rename(Sample_name = anonymized_name,
          Accession = Run)

  
  zip_list <- unzip(scale_zip, list = TRUE)
  scale_xlsx <- zip_list$Name[1]  
  temp_dir <- tempdir()
  unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
  scale_path <- file.path(temp_dir, scale_xlsx)
  
  raw <- read_xlsx(scale_path, sheet = 3,
                 col_names = FALSE)
  starts <- which(raw[[1]] %in% c("Sample ID", "Variable 1"))
  ends   <- c(starts[-1] - 1, nrow(raw))
  blocks <- map2(starts, ends, ~ raw[.x:.y, ])
  dfs <- map(blocks, function(block) {
    hdr <- block[1, ] %>% unlist() %>% as.character()
    dat <- block[-1, , drop=FALSE]
    colnames(dat) <- hdr
    dat
  })
  df1 <- dfs[[1]][, !is.na(names(dfs[[1]])), drop = FALSE]
  df2 <- as.data.frame(dfs[[2]])
  scale_s3 <- full_join(df1, df2, by = "Sample ID")
  scale_s3 <- scale_s3 %>%
    rename(Sample_name = `Sample ID`)
  scale_s3 <- scale_s3[-c(17:21), ]

  raw <- read_xlsx(scale_path, sheet = 8, col_names = FALSE)
  starts <- which(raw[[1]] %in% c("Sample ID", "Variable 1"))
  ends   <- c(starts[-1] - 1, nrow(raw))
  blocks <- map2(starts, ends, ~ raw[.x:.y, ])
  dfs <- map(blocks, function(block) {
    hdr <- block[1, ] %>% unlist() %>% as.character()
    dat <- block[-1, , drop=FALSE]
    colnames(dat) <- hdr
    dat
  })
  scale_s8 <- as.data.frame(dfs[[1]]) 
  scale_s8 <- scale_s8[-c(33:34), ]
  scale_s8 <- scale_s8 %>%
    select(
      Sample_name = `Sample ID`,
      ddPCR_Mean  = `ddPCR average copies/uL DNA`,
      ddPCR_SD    = `ddPCR SD copies/ul DNA`,
      qPCR_Mean   = `qPCR copies/ul DNA`
    )

  
  scale <- scale_s3 %>%
    full_join(scale_s8, by = "Sample_name") 

  scale <- scale %>%
    rename_with(~ case_when(
      .x == "Sample_name" ~ "Sample_name",

      str_detect(.x, "FACS cell count.*#1")  ~ "facs_rep1",
      str_detect(.x, "FACS cell count.*#2")  ~ "facs_rep2",
      str_detect(.x, "Average FACS")         ~ "facs_mean",
      str_detect(.x, "S.D.*\\.x|S.D\\.\\.x") ~ "facs_sd",

      str_detect(.x, "Ct-value replicate 1") ~ "qpcr_ct_rep1",
      str_detect(.x, "Ct-value replicate 2") ~ "qpcr_ct_rep2",
      str_detect(.x, "Average Ct-value")     ~ "qpcr_ct_mean",
      str_detect(.x, "S.D.*\\.y|S.D\\.\\.y") ~ "qpcr_sd",
      str_detect(.x, "log copies per gram")  ~ "qpcr_log",
      str_detect(.x, "Copies per gram")      ~ "qpcr_copies",

      .x == "ddPCR_Mean"                     ~ "ddpcr_mean",
      .x == "ddPCR_SD"                       ~ "ddpcr_sd",
      .x == "qPCR_Mean"                      ~ "qpcr_total_mean",

      TRUE ~ .x 
    ))

  unlink(temp_dir, recursive = TRUE)

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