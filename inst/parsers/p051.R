parse_2024_kruger_scientificreports_ddpcrhealthysubjects <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tibble", "tidyverse")
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

  library(tibble)
  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2024_kruger_scientificreports_ddpcrhealthysubjects")

  # ----- File paths -----
  counts_zip           <- file.path(local, "41598_2024_75477_MOESM2_ESM.csv.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  metadata_zip_1       <- file.path(local, "41598_2024_75477_MOESM3_ESM.csv.zip")
  metadata_zip_2       <- file.path(local, "41598_2024_75477_MOESM4_ESM.csv.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA1162476_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1162476_dada2_taxa.rds.zip")
  dnaconcentrations_zip <- file.path(local, "dnaconcentrations.csv.zip")

  subjects <- tribble(
    ~SRA_subject_ID, ~Volunteer,     ~Manuscript_subject_ID,
                7,  "Volunteer7",                     1,
                18,  "Volunteer18",                    2,
                19,  "Volunteer19",                    3,
                22,  "Volunteer22",                    4,
                36,  "Volunteer36",                    5,
                52,  "Volunteer52",                    6,
                67,  "Volunteer67",                    7,
                68,  "Volunteer68",                    8,
                88,  "Volunteer88",                    9,
                99,  "Volunteer99",                   10
  )

  # the 3 milling types
  tp_milling <- tribble(
    ~Timepoint, ~Milling,
            1,  "min",
            2,  "min",
            3,  "min",
            1,  "plus",
            2,  "plus",
            3,  "plus",
            1,  "PCRv2"
  )

  # 1) Generate every real sample (no IKA suffix here)
  base_map <- subjects %>%
    crossing(tp_milling) %>%
    mutate(
      SRA_Sample_ID = sprintf("S%s.%d.%s",
                              SRA_subject_ID, Timepoint, Milling),
      New_sample_ID = sprintf("S%s.%d.%s",
                              Manuscript_subject_ID, Timepoint, Milling)
    )

  ika_map <- base_map %>%
    filter(Timepoint == 1, Milling == "min") %>%
    mutate(
      SRA_Sample_ID = str_replace(SRA_Sample_ID, "\\.min$", ".IKA.min"),
      New_sample_ID = str_replace(New_sample_ID, "\\.min$", ".IKA.min")
    )

  # 3) Combine them
  tidy_map <- bind_rows(base_map, ika_map) %>%
    arrange(SRA_subject_ID, Milling, Timepoint)

  
  sra = read_zipped_table(metadata_zip, row.names=NULL) %>% 
      rename(Accession = Run, SRA_Sample_ID = `Sample Name`) %>%
      mutate(
        Subject = paste0(gsub(".*Volunteer(\\d+).*", "\\1", `description:_replicate`)),
        Timepoint = gsub("^S\\d+\\.(\\d+).*", "\\1", SRA_Sample_ID),
        Replicate = ifelse(grepl("Replicate", `description:_replicate`), Timepoint, 0),
        Processing = case_when(
          grepl("mill homogenized", `description:_sample_information`) ~ "mill_homogenized",
          grepl("only hammered not milled", `description:_sample_information`) ~ "hammered only",
          TRUE ~ NA_character_
        )
      )

  sra <- sra %>% 
    mutate(
      SRA_Sample_ID = case_when(

        !str_detect(SRA_Sample_ID, "min|plus|PCRv2") ~ 
          str_c(SRA_Sample_ID, ".plus"),
        # 3) otherwise leave as‐is (min and existing IKA.min / IKA.plus)
        TRUE ~ SRA_Sample_ID
      ),
      # then pull out your Milling factor in one go
      Milling = str_extract(SRA_Sample_ID, "(min|plus|PCRv2)$") %>%
                  factor(levels = c("min","plus","PCRv2"))
    )

  sra = sra %>% left_join(tidy_map, by = "SRA_Sample_ID") %>% select(-Milling.y, Timepoint.y) %>% rename(Milling = Milling.x, Timepoint = Timepoint.x, SampleID = New_sample_ID)

  # ---- counts ----
  if (!file.exists(counts_zip)) stop("Counts file not found: ", counts_zip)
  dataset <- read_zipped_table(counts_zip, row.names = NULL) %>% 
    mutate(Frequency = gsub("Day ", "", Frequency)) %>%
    mutate(Milling = ifelse(Milling == "Min", "PCRv2", tolower(Milling))) %>%
    mutate(SampleID = sprintf("S%s.%d.%s", Subject, Timepoint, Milling))
  colnames(dataset) <- ifelse(
    grepl("^[kpcofg]__(?!.*__).*", colnames(dataset), perl = TRUE),
    gsub("([kpcofg])__", "\\1_", colnames(dataset)),
    colnames(dataset)
  )

  # ---- metadata ----
  metadata_cols <- c("SampleID", "Subject", "Timepoint", "Milling", "Frequency", "StoolsperDay",
                     "BristolStoolScale_highest", "WaterContent_perc", "pH",
                     "Calprotectin_ugperg", "MPO_ngperml", "PhylogeneticDiversity", "Chao1", "inverse_simpson", "gini_simpson", "shannon", "fisher")
  metadata <- dataset[, metadata_cols] 
  
  metadata1 <- read_zipped_table(metadata_zip_1, row.names = NULL) %>%
    mutate(Replicate = gsub("Replicate ", "Replicate_", ID))
  metadata2 <- read_zipped_table(metadata_zip_2, row.names = NULL) %>%
    pivot_longer(cols = starts_with("Replicate"), names_to = "Replicate", values_to = "Value") %>%
    pivot_wider(names_from = ID, values_from = Value)

  merged_metadata <- metadata1 %>%
    full_join(metadata2, by = "Replicate") 

  metadata <- metadata %>%
  mutate(
    SampleID = factor(SampleID),
    Subject                     = factor(Subject),
    Timepoint                   = factor(Timepoint)
  ) %>%
  
  mutate(
    Frequency     = as.integer(Frequency),
    StoolsperDay  = as.integer(StoolsperDay)
  ) %>%
  
  mutate(
    Milling = factor(Milling)
  ) %>%
  
  mutate(
    BristolStoolScale_highest = factor(
      as.integer(BristolStoolScale_highest),
      levels   = 1:6,
      ordered  = TRUE
    ),

    ## 4. pH ─ numeric
    pH = as.numeric(pH),
    Chao1                     = as.integer(Chao1),
    across(
      c(WaterContent_perc, pH,
        Calprotectin_ugperg, MPO_ngperml,
        PhylogeneticDiversity,
        inverse_simpson, gini_simpson,
        shannon, fisher),
      as.numeric
    )
  ) %>% 
  mutate(Timepoint = factor(Timepoint,
                            levels   = sort(unique(Timepoint)),
                            ordered  = TRUE)) 

  metadata = full_join(metadata, sra, by = "SampleID")

  dnaconcentrations = read_zipped_table(dnaconcentrations_zip, row.names = NULL) 
  metadata = full_join(metadata, dnaconcentrations, by = "SRA_Sample_ID") %>% rename(Milling = Milling.y, Subject = Subject.y, Timepoint = Timepoint.y.y) %>% select(-c(Milling.x, Subject.x, Timepoint.x))
    
  # ---- scale ----
  scale_cols <- c("SampleID",
                  "Mean Fungi copies_per mg total weight", 
                  "Mean Fungi copies_per mg dry weight", 
                  "Mean bacteria copies_per mg total weight", 
                  "Mean bacteria copies_per mg dry weight")
  scale    <- dataset[, scale_cols] 
  scale <- scale %>% mutate(log2_Mean_Fungi_copies_per_mg_total_weight = ifelse(`Mean Fungi copies_per mg total weight` > 0, log2(`Mean Fungi copies_per mg total weight`),NA),
                         log2_Mean_Fungi_copies_per_mg_dry_weight = ifelse(`Mean Fungi copies_per mg dry weight` > 0, log2(`Mean Fungi copies_per mg dry weight`),NA),
                         log2_Mean_bacteria_copies_per_mg_total_weight = ifelse(`Mean bacteria copies_per mg total weight` > 0, log2(`Mean bacteria copies_per mg total weight`),NA),
                         log2_Mean_bacteria_copies_per_mg_dry_weight = ifelse(`Mean bacteria copies_per mg dry weight` > 0, log2(`Mean bacteria copies_per mg dry weight`),NA),
                         log10_Mean_Fungi_copies_per_mg_total_weight = ifelse(`Mean Fungi copies_per mg total weight` > 0, log10(`Mean Fungi copies_per mg total weight`),NA),
                         log10_Mean_Fungi_copies_per_mg_dry_weight = ifelse(`Mean Fungi copies_per mg dry weight` > 0, log10(`Mean Fungi copies_per mg dry weight`),NA),
                         log10_Mean_bacteria_copies_per_mg_total_weight = ifelse(`Mean bacteria copies_per mg total weight` > 0, log10(`Mean bacteria copies_per mg total weight`),NA),
                         log10_Mean_bacteria_copies_per_mg_dry_weight = ifelse(`Mean bacteria copies_per mg dry weight` > 0, log10(`Mean bacteria copies_per mg dry weight`),NA)) 

  scale <- scale %>% left_join(sra %>% select(SampleID), by = "SampleID") %>%
  filter(!if_all(.cols = -SampleID, .fns = is.na))
  
  scale_cols <- c("SampleID", "TotalSCFAs_umolpermg", "TotalBCFAs_umolpermg" )
  metabolomicsscale <- dataset %>% 
    select(scale_cols) %>%
    mutate(log2_TotalSCFAs_umolpermg = ifelse(`TotalSCFAs_umolpermg` > 0, log2(`TotalSCFAs_umolpermg`),NA),
           log2_TotalBCFAs_umolpermg = ifelse(`TotalBCFAs_umolpermg` > 0, log2(`TotalBCFAs_umolpermg`),NA),
           log10_TotalSCFAs_umolpermg = ifelse(`TotalSCFAs_umolpermg` > 0, log10(`TotalSCFAs_umolpermg`),NA),
           log10_TotalBCFAs_umolpermg = ifelse(`TotalBCFAs_umolpermg` > 0, log10(`TotalBCFAs_umolpermg`),NA))

  counts_original <- dataset %>% 
    select(-c(3:31)) %>%
    column_to_rownames("SampleID") %>%
    select(-c("Subject", "Timepoint"))
    
  counts_metabolomics <- dataset %>% 
    select(SampleID, c(13:23)) %>%
    column_to_rownames("SampleID")

  # ---- tax ----
  tax <- tibble(
    taxonomy = colnames(counts_original)[!(colnames(counts_original) %in% c(metadata_cols, scale_cols))]
  )

  tax_metabolomics <- tibble(
    metabolites = colnames(counts_metabolomics)[!(colnames(counts_metabolomics) %in% c(metadata_cols, scale_cols))]
  )

  parse_taxonomy_column <- function(name) {
    # helper: capitalize first letter
    cap1 <- function(x) {
      if (is.na(x) || x == "") return(NA_character_)
      paste0(toupper(substr(x,1,1)), tolower(substr(x,2,nchar(x))))
    }

    # prepare output
    ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","ogtaxonomy")
    out   <- setNames(as.list(rep(NA_character_, length(ranks))), ranks)
    out$ogtaxonomy <- name

    if (grepl("^k__", name)) {
      parts <- strsplit(name, "__", fixed = TRUE)[[1]]
      # parts[1] is like "k", parts[2] like "Bacteria_p", parts[3] like "Actinobacteriota_c", …
      for (i in seq(2, length(parts))) {
        # prefix code is last char of parts[i-1], e.g. "k", "p", "c", …
        prefix <- sub(".*_", "", parts[i-1])
        # value is parts[i] with any trailing "_<code>" removed
        value  <- gsub("_[kpcofg]$", "", parts[i])
        value  <- gsub("_", " ", value)   # optional: turn underscores to spaces
        value  <- cap1(value)
        switch(prefix,
          k = out$Kingdom  <- value,
          p = out$Phylum   <- value,
          c = out$Class    <- value,
          o = out$Order    <- value,
          f = out$Family   <- value,
          g = out$Genus    <- value,
          NULL
        )
      }
    } else if (grepl("^g_{1,2}", name)) {
      # short genus label
      genus <- sub("^g_{1,2}", "", name)
      genus <- gsub("_", " ", genus)
      out$Genus <- cap1(genus)
    }

    out
  }

  # Apply to your counts_original columns:
  taxonomy_df <- do.call(
    rbind,
    lapply(colnames(counts_original), parse_taxonomy_column)
  )
  rownames(taxonomy_df) <- NULL
  taxonomy_df <- as.data.frame(taxonomy_df, stringsAsFactors = FALSE)
  taxonomy_df = make_taxa_label(taxonomy_df)
  rownames(taxonomy_df) <- taxonomy_df$ogtaxonomy

  if (!raw) {
    aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="SampleID", align = align, study_name=basename(local))
    counts_original = aligned$counts_original
    matched_taxa <- taxonomy_df$Taxa[match(colnames(counts_original), rownames(taxonomy_df))]
    colnames(counts_original) <- matched_taxa
    counts_original <- collapse_duplicate_columns_exact(counts_original)
    original_names <- colnames(counts_original)
    counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
    counts_original = fill_na_zero_numeric(counts_original)

    aligned = rename_and_align(counts_original = counts_metabolomics, metadata=metadata, scale=scale, by_col="SampleID", align = align, study_name=basename(local))
    counts_metabolomics = aligned$counts_original
    original_names <- colnames(counts_metabolomics)
    counts_metabolomics <- as.data.frame(lapply(counts_metabolomics, as.numeric), row.names = rownames(counts_metabolomics), col.names = original_names, check.names = FALSE)
  }

  # ---- proportions ----
  proportions <- sweep(counts_original, 1, rowSums(counts_original), '/')
  proportions_metabolomics <- sweep(counts_metabolomics, 1, rowSums(counts_metabolomics), '/')

  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA 
  tax_reprocessed2 <- NA
  counts_reprocessed2 <- NA
  proportions_reprocessed2 <- NA

  # ----- Reprocessed counts from RDS ZIP -----
  if (all(file.exists(c(repro_counts_rds_zip, repro_tax_zip)))) {
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
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="SampleID", align = align, study_name=basename(local))
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
      counts_original = fill_na_zero_numeric(counts_original)
      counts_metabolomics = fill_na_zero_numeric(counts_metabolomics)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
      proportions = fill_na_zero_numeric(proportions)
      proportions_metabolomics = fill_na_zero_numeric(proportions_metabolomics)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
  }

  return(list(
    counts      = list(original = list(metabolomics = counts_metabolomics, amplicon = counts_original), reprocessed = list(metabolomics = NA, amplicon = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2))),
    proportions = list(original = list(metabolomics = proportions_metabolomics, amplicon = proportions), reprocessed = list(metabolomics = NA, amplicon = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2))),
    tax         = list(original = list(metabolomics = tax_metabolomics, amplicon = tax), reprocessed = list(metabolomics = NA, amplicon = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2))),
    scale       = scale,
    metadata    = list(replicates = merged_metadata, metadata = metadata)
  ))
}
