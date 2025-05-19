parse_2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  if (!is.logical(raw) || length(raw) != 1) {
    stop("`raw` must be a single logical value (TRUE or FALSE)")
  }
  if (!is.logical(align) || length(align) != 1) {
    stop("`align` must be a single logical value (TRUE or FALSE)")
  }
  library(tidyverse)
  library(readxl)

  counts <- NA
  proportions <- NA
  tax <- NA
  counts_reprocessed2 <- NA
  tax_reprocessed2 <- NA
  proportions_reprocessed2 <- NA
  counts_reprocessed <- NA
  tax_reprocessed <- NA
  proportions_reprocessed <- NA
  

  # ----- Local base directory -----
  local <- file.path("2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy")

  # ----- File paths -----
  metadata_zip         <- file.path(local, "SraRunTable.csv.zip")
  scale_zip            <- file.path(local, "Table 1.XLSX.zip")
  repro_counts_rds_zip <- file.path(local, "ERP108719_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "ERP108719_dada2_taxa.rds.zip")


  # ----- Metadata -----
  metadata <- read_zipped_table(metadata_zip, row.names = NULL) %>% as.data.frame()
  metadata <- metadata %>% rename(Sample_name_existing = Sample_name, Sample_name = anonymized_name, Accession = Run)
  rownames(metadata) <- metadata$Sample_name

  # ----- Scale ----- # SHEET 9 has mock vs sample A,B,C for ddPCR replicates vs qPCR replicates at different dilution factors. We have not done that yet.
  zip_list <- unzip(scale_zip, list = TRUE)
  scale_xlsx <- zip_list$Name[1]  
  temp_dir <- tempdir("oro")
  dir.create(temp_dir)
  unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
  scale_path <- file.path(temp_dir, scale_xlsx)
  
  rawd <- read_xlsx(scale_path, sheet = 3,
                 col_names = FALSE)
  starts <- which(rawd[[1]] %in% c("Sample ID", "Variable 1"))
  ends   <- c(starts[-1] - 1, nrow(rawd))
  blocks <- map2(starts, ends, ~ rawd[.x:.y, ])
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

  rawd <- read_xlsx(scale_path, sheet = 8, col_names = FALSE)
  starts <- which(rawd[[1]] %in% c("Sample ID", "Variable 1"))
  ends   <- c(starts[-1] - 1, nrow(rawd))
  blocks <- map2(starts, ends, ~ rawd[.x:.y, ])
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
      ddPCR_copies_ul_dna_rep1 = `ddPCR copies per uL DNA replicate 1`,
      ddPCR_copies_ul_dna_rep2 = `ddPCR copies per uL DNA replicate 2`,
      qPCR_copies_ul_dna = `qPCR copies/ul DNA`
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

      str_detect(.x, "ddpcr_copies_ul_dna_rep1") ~ "ddpcr_copies_ul_dna_rep1",
      str_detect(.x, "ddpcr_copies_ul_dna_rep2") ~ "ddpcr_copies_ul_dna_rep2",
      str_detect(.x, "qpcr_copies_ul_dna") ~ "qpcr_copies_ul_dna",

      .x == "ddPCR_Mean"                     ~ "ddpcr_mean",
      .x == "ddPCR_SD"                       ~ "ddpcr_sd",
      .x == "qPCR_Mean"                      ~ "qpcr_mean",

      TRUE ~ .x 
    )) 
  scale$Sample_name <- gsub("-", ".", scale$Sample_name)
  scale<- scale %>%
  pivot_longer(
    cols        = matches("_rep[12]$"),               
    names_to    = c(".value", "Replicate"),            
    names_pattern= "(.*)_rep(\\d)$"                     
  ) %>% select(Sample_name, Replicate, facs, qpcr_ct, qpcr_log, qpcr_copies, ddPCR_copies_ul_dna, qPCR_copies_ul_dna) %>%
    mutate(across(-Sample_name, as.numeric)) %>%
    mutate(log10_FC = ifelse(facs > 0, log10(facs), NA)) %>%
    mutate(log2_FC = ifelse(facs > 0, log2(facs), NA)) %>%
    mutate(log2_qpcr_ct = ifelse(qpcr_ct > 0, log2(qpcr_ct), NA)) %>%
    mutate(log10_qpcr_ct = ifelse(qpcr_ct > 0, log10(qpcr_ct), NA)) %>%
    mutate(log2_ddpcr_copies_ul_dna = ifelse(ddPCR_copies_ul_dna > 0, log2(ddPCR_copies_ul_dna), NA)) %>%
    mutate(log10_ddpcr_copies_ul_dna = ifelse(ddPCR_copies_ul_dna > 0, log10(ddPCR_copies_ul_dna), NA)) %>%
    mutate(Sample_name = ifelse(str_detect(Sample_name, "\\.PMA$"),
                               str_replace(Sample_name, "\\.PMA$", paste0(".", Replicate, ".PMA")),
                               paste0(Sample_name, ".", Replicate))) %>%
    # Duplicate rows for sample 26 and modify their names
    bind_rows(
      filter(., str_detect(Sample_name, "^26\\.[12]\\.PMA$")) %>%
        mutate(Sample_name = str_replace(Sample_name, "\\.PMA$", ".duplo.PMA"))
    )

  cleanup_tempfiles(temp_dir)

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
    temp_dir <- tempdir("repro")
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
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample_name", align = align, study_name = basename(local))
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
  }

  if (!raw) {
    counts = fill_na_zero_numeric(counts)
    counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
    proportions = fill_na_zero_numeric(proportions)
    proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
    proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
  }

  return(list(
    counts = list(
                original = counts,
                reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
    ),
    proportions = list(
                original = proportions,
                reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
    ),
    tax = list(
                original = tax,
                reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
    ),
    scale = scale,
    metadata = metadata
  ))
}