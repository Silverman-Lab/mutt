parse_2024_tunsakul_peerj_aerobicvsanaerobicinhealthyvsobesity <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse")
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
  
  # ----- Local base directory -----
  local <- file.path("2024_tunsakul_peerj_aerobicvsanaerobicinhealthyvsobesity")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "PRJNA1020208_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1020208_dada2_taxa.rds.zip")
  metadata_zip         <- file.path(local, "SraRunTable (38).csv.zip")
  scale_zip            <- file.path(local, "scale.csv.zip") 

  # --- original counts, proportions, tax ---
  counts <- NA
  proportions <- NA
  tax <- NA
  
  # --- metadata and scale ----
  metadata <- read_zipped_table(metadata_zip, row.names = NULL) %>%
    as.data.frame() %>%
    rename(Sample = `Sample Name`, Accession = Run) %>%
    mutate(
      Sample = as.character(Sample),
      
      # Remove any existing 'a', 'aa', or 'anan' suffix
      Sample = gsub("(a|aa|anan|an)$", "", Sample),
      
      # Add ID prefix if it's numeric
      Sample = ifelse(grepl("^\\d+$", Sample), paste0("ID", Sample), Sample),
      
      # Append new suffix based on collection method
      Sample = paste0(Sample, if_else(samp_collect_device == "aerobic collection", "a", "an"))
    ) %>%
    rename(environment = samp_collect_device)



  suppdoc      <- tribble(
    ~Sample, ~Type, ~Raw_Reads, ~Quality_Reads,
    "ID1a",   "a", 118505,  94439,
    "ID2a",   "a", 80082,   60357,
    "ID3a",   "a", 71065,   48416,
    "ID4a",   "a", 104163,  78155,
    "ID5a",   "a", 81255,   58685,
    "ID6a",   "a", 97153,   67858,
    "ID7a",   "a", 89970,   59716,
    "ID8a",   "a", 79423,   60230,
    "ID9a",   "a", 84451,   54372,
    "ID10a",  "a", 91965,   65186,
    "ID11a",  "a", 84075,   51981,
    "ID12a",  "a", 108682,  71907,
    "ID13a",  "a", 88285,   64335,
    "ID14a",  "a", 26181,   16724,
    "ID15a",  "a", 20825,   14460,
    "ID16a",  "a", 120845,  87511,
    "ID17a",  "a", 84156,   54143,
    "ID18a",  "a", 22177,   13904,
    "ID19a",  "a", 36802,   23328,
    "ID20a",  "a", 27583,   16628,
    "ID1an",  "an", 22743,  16527,
    "ID2an",  "an", 25254,  15573,
    "ID3an",  "an", 11463,   7137,
    "ID4an",  "an", 32891,  21548,
    "ID5an",  "an", 34344,  24527,
    "ID6an",  "an", 50952,  33129,
    "ID7an",  "an", 119967, 78021,
    "ID8an",  "an", 60443,  42347,
    "ID9an",  "an", 67837,  40939,
    "ID10an", "an", 43606,  30408,
    "ID11an", "an", 48478,  31043,
    "ID12an", "an", 48462,  29440,
    "ID13an", "an", 21089,  12679,
    "ID14an", "an", 21849,  12940,
    "ID15an", "an", 29929,  18629,
    "ID16an", "an", 38978,  29505,
    "ID17an", "an", 64451,  46372,
    "ID18an", "an", 63622,  44519,
    "ID19an", "an", 22446,  13066,
    "ID20an", "an", 19512,  12833
  )      

  suppdoc <- suppdoc %>% select(-Type)

  metadata <- metadata %>%
    left_join(suppdoc, by = "Sample")

  scale        <- read_zipped_table(scale_zip, row.names=NULL) %>% as.data.frame() %>% 
                  rename(`qPCR (copies/g)` = `Approximate Mean Load`) %>%
                  mutate(log2_qPCR = ifelse(`qPCR (copies/g)` > 0,log2(`qPCR (copies/g)`), NA)) %>%
                  mutate(log10_qPCR = ifelse(`qPCR (copies/g)` > 0, log10(`qPCR (copies/g)`), NA)) %>%
                  mutate(Sample = paste0("ID",Sample)) %>%
                  select(-Condition)

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
    temp_dir <- tempfile("repro")
    dir.create(temp_dir)
    unzipped <- unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
    counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(counts_file))

    # ----- Taxonomy reprocessed -----
    unzipped <- unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))

    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name = basename(local))
      counts_reprocessed = aligned$reprocessed
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
      original_names <- colnames(counts_reprocessed)
      counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
    }
    
    # proportions reprocessed
    proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
    cleanup_tempfiles(temp_dir)
  }

  if (!raw) {
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }
  
  return(list(
    counts = list(
      original = NA,
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original = NA,
      reprocessed = proportions_reprocessed
    ),
    tax = list(
      original = NA,
      reprocessed = tax_reprocessed
    ),
    scale = scale,
    metadata = metadata
  ))
}
