parse_2021_marotz_mSystems_oral_mouthwash <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse")
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

  # ----- Local base directory -----
  local <- file.path("2021_marotz_mSystems_oral_mouthwash")

  # ----- File paths -----
  counts_zip    <- file.path(local, "2021_marotz_mSystems_oral_mouthwash.RDS.zip")
  metadata_zip  <- file.path(local, "T3_SRS_metadata_ms.txt.zip")
  sra_zip       <- file.path(local, "SraRunTable (40).csv.zip")
  sra_zip2       <- file.path(local, "SraRunTable (47).csv.zip")



  repro_counts_zips <- c(
    file.path(local, "ERP111447_dada2_counts.rds.zip"),
    file.path(local, "ERP117149_dada2_counts.rds.zip")
  )

  repro_tax_zips <- c(
    file.path(local, "ERP111447_dada2_taxa.rds.zip"),
    file.path(local, "ERP117149_dada2_taxa.rds.zip")
  )

  read_zipped_tabled <- function(zip_path, sep = ",", header = TRUE,
                                row.names = 1, check.names = FALSE,
                                skip = 0) {
    inner_file <- unzip(zip_path, list = TRUE)$Name[1]
    con        <- unz(zip_path, inner_file)
    read.table(con,
              sep          = sep,
              header       = header,
              row.names    = row.names,
              check.names  = check.names,
              stringsAsFactors = FALSE,
              skip         = skip,
              fill         = TRUE,
              comment.char = "",
              quote        = "\"")
  }


  # ----- Metadata and Scale -----
  metadata <- read_zipped_table(metadata_zip, sep=",", row.names = NULL)
  sra1 = read_zipped_table(sra_zip, row.names = NULL)
  sra1 <- sra1 %>% rename(Accession = Run, saliva_sample_ID = saliva_sample_id, 
         FC_avg_cells_5_min = fc_avg_cells_5_min, FC_avg_cells_per_ul = fc_avg_cells_per_ul,
         FC_cells_per_ul_r1 = fc_cells_per_ul_r1, FC_cells_per_ul_r2 = fc_cells_per_ul_r2) 
  sra1$SampleID <- paste0(sra1$saliva_sample_ID, ".", sra1$timepoint, ".", sra1$processing)
  bad_values <- c(".NA.", paste0(".", 1:9, "."))
  sra1$SampleID[sra1$SampleID %in% bad_values] <- NA
  sra1 <- sra1 %>%
    mutate(
      SampleID = if_else(
        `Assay Type` == "WGS" & is.na(SampleID),
        {
          wd_clean <- str_remove(`well_description (exp)`, "\\.$")
          parts     <- str_match(wd_clean, "^(.*?)\\.ly(.*)$")
          prefix    <- parts[,2]
          suffix    <- parts[,3]
          paste0(prefix, ".", timepoint, ".", suffix)
        },
        SampleID
      )
    ) %>%
  mutate(
    SampleID = if_else(
      is.na(SampleID),
      str_replace(`anonymized_name`, "\\.ly", "."),
      SampleID
    )
  ) %>%
    mutate(
      SampleID = if_else(
        SampleID %in% c("NA.NA.NA", "") | is.na(SampleID),
        anonymized_name,
        SampleID
      )
    ) %>% mutate(SampleID = if_else(
                          SampleID %in% c("") | is.na(SampleID),
                          host_subject_id, SampleID)) %>%
  mutate(
    SampleID = if_else(
      str_length(SampleID) == 1,
      {
        parts <- str_match(`orig_name (exp)`, "^(.*?)\\.(.*)$")
        prefix <- parts[,2]
        suffix <- parts[,3]
        paste0(prefix, ".", timepoint, ".", suffix)
      },
      SampleID
    )
  ) %>%
    mutate(
      SampleID = if_else(
        str_detect(SampleID, "BLANK") &
          !is.na(`orig_name (exp)`) &
          `orig_name (exp)` != "",
        `orig_name (exp)`,
        SampleID
      )
    ) 
  sra2 = read_zipped_tabled(sra_zip2, row.names = NULL) %>% 
  rename(Accession = Run, saliva_sample_ID = saliva_sample_id, 
         FC_avg_cells_5_min = fc_avg_cells_5_min, FC_avg_cells_per_ul = fc_avg_cells_per_ul,
         FC_cells_per_ul_r1 = fc_cells_per_ul_r1, FC_cells_per_ul_r2 = fc_cells_per_ul_r2) %>%
  mutate(
    SampleID = paste0(saliva_sample_ID, ".", processing),
    SampleID = if_else(
      SampleID == ".",
      anonymized_name,
      SampleID
    )
  )

  common_cols <- intersect(names(sra1), names(sra2))
  types1 <- sapply(sra1[common_cols], function(x) class(x)[1])
  types2 <- sapply(sra2[common_cols], function(x) class(x)[1])
  mismatch <- common_cols[types1 != types2]
  sra1 <- sra1 %>% mutate(across(all_of(mismatch), as.character))
  sra2 <- sra2 %>% mutate(across(all_of(mismatch), as.character))

  merge_same_cols <- function(df) {
    cols_x <- grep("\\.x$", names(df), value = TRUE)
    for (col_x in cols_x) {
      base   <- sub("\\.x$", "", col_x)
      col_y  <- paste0(base, ".y")
      if (col_y %in% names(df)) {
        x_vals <- df[[col_x]]
        y_vals <- df[[col_y]]
        conflict <- !is.na(x_vals) & !is.na(y_vals) & x_vals != y_vals
        if (any(conflict)) {
          warning(sprintf(
            "Column '%s': %d conflicting rows (first at row %d). Keeping the .x version there.",
            base, sum(conflict), which(conflict)[1]
          ))
        }
        df[[base]] <- coalesce(x_vals, y_vals)
        df[[col_x]] <- NULL
        df[[col_y]] <- NULL
      }
    }
    df
  }

  sra <- full_join(sra1,sra2,by      = "Accession",suffix  = c(".x", ".y")) %>% merge_same_cols()
  sra$SampleID <- ifelse(
    sra$`Assay Type` == "WGS",
    paste0(sra$SampleID, ".", sra$`Assay Type`),
    sra$SampleID
  )
  sra <- type.convert(sra, as.is = TRUE)
  sra = remove_empty_columns(sra)
  sra$SampleID <- paste0(sra$SampleID, ".", sra$`run_date (exp)`)

  metadata <- metadata %>% rename(Sample = SampleID, Participant_ID = participant_id)
  metadata$SampleID <- paste0(metadata$saliva_sample_ID, ".", metadata$processing)
  metadata <- type.convert(metadata, as.is = TRUE)
  #metadata <- full_join(sra, metadata, by = c("Participant_ID", "saliva_weight_g", "FC_cells_per_ul_r1", "FC_cells_per_ul_r2", "FC_avg_cells_per_ul", "FC_avg_cells_5_min"))
  metadata = remove_empty_columns(metadata)

  scale_sra <- sra %>% select(SampleID, FC_cells_per_ul_r1, FC_cells_per_ul_r2, FC_avg_cells_per_ul, 
                              FC_avg_cells_5_min, qpcr_median_16s_copies_per_2ul_dna, all_flow_cells_5min_avg,
                              all_flow_cellsperul_avg, all_qpcr_cells_5min_avg, all_qpcr_cellsperul_avg,
                              live_flow_cells_5min_avg, live_flow_cellsperul_avg, live_qpcr_cells_5min_avg,
                              live_qpcr_cellsperul_avg) %>% 
           mutate(FC_sd_cells_per_ul = sqrt((FC_cells_per_ul_r1 - FC_avg_cells_per_ul)^2 + (FC_cells_per_ul_r2 - FC_avg_cells_per_ul)^2) / 2)

  scale_sra = scale_sra %>% 
    mutate(log2_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log2(FC_avg_cells_per_ul), NA)) %>%
    mutate(log10_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log10(FC_avg_cells_per_ul), NA)) %>%
    mutate(log2_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log2(FC_avg_cells_5_min), NA)) %>%
    mutate(log10_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log10(FC_avg_cells_5_min), NA)) %>%
    mutate(log2_FC_sd_cells_per_ul = ifelse(FC_sd_cells_per_ul > 0, log2(FC_sd_cells_per_ul), NA)) %>%
    mutate(log10_FC_sd_cells_per_ul = ifelse(FC_sd_cells_per_ul > 0, log10(FC_sd_cells_per_ul), NA)) %>%
    mutate(log2_qpcr_median_16s_copies_per_2ul_dna = ifelse(qpcr_median_16s_copies_per_2ul_dna > 0, log2(qpcr_median_16s_copies_per_2ul_dna), NA)) %>%
    mutate(log10_qpcr_median_16s_copies_per_2ul_dna = ifelse(qpcr_median_16s_copies_per_2ul_dna > 0, log10(qpcr_median_16s_copies_per_2ul_dna), NA)) %>%
    mutate(log2_all_flow_cells_5min_avg = ifelse(all_flow_cells_5min_avg > 0, log2(all_flow_cells_5min_avg), NA)) %>%
    mutate(log10_all_flow_cells_5min_avg = ifelse(all_flow_cells_5min_avg > 0, log10(all_flow_cells_5min_avg), NA)) %>%
    mutate(log2_all_flow_cellsperul_avg = ifelse(all_flow_cellsperul_avg > 0, log2(all_flow_cellsperul_avg), NA)) %>%
    mutate(log10_all_flow_cellsperul_avg = ifelse(all_flow_cellsperul_avg > 0, log10(all_flow_cellsperul_avg), NA)) %>%
    mutate(log2_all_qpcr_cells_5min_avg = ifelse(all_qpcr_cells_5min_avg > 0, log2(all_qpcr_cells_5min_avg), NA)) %>%
    mutate(log10_all_qpcr_cells_5min_avg = ifelse(all_qpcr_cells_5min_avg > 0, log10(all_qpcr_cells_5min_avg), NA)) %>%
    mutate(log2_all_qpcr_cellsperul_avg = ifelse(all_qpcr_cellsperul_avg > 0, log2(all_qpcr_cellsperul_avg), NA)) %>%
    mutate(log10_all_qpcr_cellsperul_avg = ifelse(all_qpcr_cellsperul_avg > 0, log10(all_qpcr_cellsperul_avg), NA)) %>%
    mutate(log2_live_flow_cells_5min_avg = ifelse(live_flow_cells_5min_avg > 0, log2(live_flow_cells_5min_avg), NA)) %>%
    mutate(log10_live_flow_cells_5min_avg = ifelse(live_flow_cells_5min_avg > 0, log10(live_flow_cells_5min_avg), NA)) %>%
    mutate(log2_live_flow_cellsperul_avg = ifelse(live_flow_cellsperul_avg > 0, log2(live_flow_cellsperul_avg), NA)) %>%
    mutate(log10_live_flow_cellsperul_avg = ifelse(live_flow_cellsperul_avg > 0, log10(live_flow_cellsperul_avg), NA)) 

  scale_metadata <- metadata %>% select(SampleID, FC_cells_per_ul_r1, FC_cells_per_ul_r2, FC_avg_cells_per_ul, 
                              FC_avg_cells_5_min) %>% 
           mutate(FC_sd_cells_per_ul = sqrt((FC_cells_per_ul_r1 - FC_avg_cells_per_ul)^2 + (FC_cells_per_ul_r2 - FC_avg_cells_per_ul)^2) / 2)

  scale_metadata = scale_metadata %>% 
    mutate(log2_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log2(FC_avg_cells_per_ul), NA)) %>%
    mutate(log10_FC_avg_cells_per_ul = ifelse(FC_avg_cells_per_ul > 0, log10(FC_avg_cells_per_ul), NA)) %>%
    mutate(log2_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log2(FC_avg_cells_5_min), NA)) %>%
    mutate(log10_FC_avg_cells_5_min = ifelse(FC_avg_cells_5_min > 0, log10(FC_avg_cells_5_min), NA)) %>%
    mutate(log2_FC_sd_cells_per_ul = ifelse(FC_sd_cells_per_ul > 0, log2(FC_sd_cells_per_ul), NA)) %>%
    mutate(log10_FC_sd_cells_per_ul = ifelse(FC_sd_cells_per_ul > 0, log10(FC_sd_cells_per_ul), NA))

  # ----- Original Counts and Taxonomy -----
  orig_rds_file <- unzip(counts_zip, list = TRUE)$Name[1]
  temp_dir <- tempfile("repro")
  dir.create(temp_dir)
  unzip(counts_zip, files = orig_rds_file, exdir = temp_dir, overwrite = TRUE)
  counts_original <- readRDS(file.path(temp_dir, orig_rds_file)) %>% t() %>% as.data.frame()
  counts_original <- collapse_duplicate_columns_exact(counts_original)
  cleanup_tempfiles(temp_dir)
  counts_original <- counts_original %>%
    rownames_to_column("SampleID") %>%
    mutate(SampleID = if_else(
      str_detect(SampleID, "raw"),
      SampleID,
      paste0(SampleID, ".PMA")
    ))
  counts_original <- counts_original %>% column_to_rownames("SampleID")

  raw_tax <- data.frame(Taxa = colnames(counts_original))
  tax_original <- raw_tax %>%
    mutate(taxa = str_trim(Taxa)) %>%
    separate(
      Taxa,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = "\\s*;\\s*",
      extra = "drop",
      fill = "right"
    )
  tax_original$species = tax_original$Species
  tax_original$Species = NULL
  tax_original = make_taxa_label(tax_original)
  tax_original$Species = tax_original$species
  tax_original$species = NULL
  rownames(tax_original) <- tax_original$taxa

  if (!raw) {
    aligned = rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale_metadata, by_col = "SampleID", align = align, study_name = basename(local))
    counts_original = aligned$counts_original
    matched_taxa <- tax_original$Taxa[match(colnames(counts_original), rownames(tax_original))]
    colnames(counts_original) <- matched_taxa
    counts_original = collapse_duplicate_columns_exact(counts_original)
    original_names <- colnames(counts_original)
    counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
  }
  proportions_original <- sweep(counts_original, 1, rowSums(counts_original), "/")


  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA
  counts_reprocessed2 <- NA
  proportions_reprocessed2 <- NA
  tax_reprocessed2 <- NA

  if (all(file.exists(repro_counts_zips), file.exists(repro_tax_zips))) {
    # Process multiple zipped RDS files
    counts_reprocessed_list <- list()
    proportions_reprocessed_list <- list()
    tax_reprocessed_list <- list()
    counts_reprocessed2_list <- list()
    proportions_reprocessed2_list <- list()
    tax_reprocessed2_list <- list()

    for (i in seq_along(repro_counts_zips)) {
        # Unzip and read counts
        temp_rds <- tempfile("repro")
        dir.create(temp_rds)
        unzipped = unzip(repro_counts_zips[i], exdir = temp_rds, overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop(paste("No *_counts.rds file found for index", i))
        counts <- as.data.frame(readRDS(counts_file))

        # ----- rdp16 -----
        if (!file.exists(file.path(local,"rdp16classified.csv.zip"))) {
          if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
              required_pkgs <- c("dada2", "Biostrings")
              missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
              if (length(missing_pkgs) > 0) {
                stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                    ". Please install them before running this function.")
              }
              seqs <- Biostrings::DNAStringSet(colnames(counts))
              rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
              tax_reprocessed2 = make_taxa_label(rdpclassified) 
            } else {
              stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
          }
          
          } else {
            tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
        }

        # Unzip and read taxonomy
        unzipped = unzip(repro_tax_zips[i], exdir = temp_rds, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop(paste("No *_taxa.rds file found for index", i))
        tax <- as.data.frame(readRDS(tax_file))
        tax <- make_taxa_label(tax)

        if (!raw) {
          aligned = rename_and_align(counts_reprocessed = counts, metadata = sra, scale = scale_sra, by_col = "SampleID", align = align, study_name = basename(local))
          counts = aligned$reprocessed
          if (nrow(counts) > 0) {
            matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
            colnames(counts) <- matched_taxa
            counts = collapse_duplicate_columns_exact(counts)
            original_names <- colnames(counts)
            counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
          }
          counts2 = aligned$reprocessed
          if (nrow(counts2) > 0) {
            matched_taxa <- tax_reprocessed2$Taxa[match(colnames(counts2), rownames(tax_reprocessed2))]
            colnames(counts2) <- matched_taxa
            counts2 = collapse_duplicate_columns_exact(counts2)
            original_names <- colnames(counts2)
            counts2 <- as.data.frame(lapply(counts2, as.numeric), row.names = rownames(counts2), col.names = original_names, check.names = FALSE)
          }
        }

        if (nrow(counts) > 0) {
          # proportions
          proportions <- sweep(counts, 1, rowSums(counts), '/')
          proportions2 <- sweep(counts2, 1, rowSums(counts2), '/')
          # Label with study name based on zip filename prefix
          study_id <- sub("_.*$", "", basename(tools::file_path_sans_ext(repro_counts_zips[i])))
          counts$Study <- study_id
          proportions$Study <- study_id
          tax$Study <- study_id
          counts2$Study <- study_id
          proportions2$Study <- study_id
          tax_reprocessed2$Study <- study_id

          counts_reprocessed_list[[i]] <- counts
          proportions_reprocessed_list[[i]] <- proportions
          tax_reprocessed_list[[i]] <- tax
          counts_reprocessed2_list[[i]] <- counts2
          proportions_reprocessed2_list[[i]] <- proportions2
          tax_reprocessed2_list[[i]] <- tax_reprocessed2
        }

        cleanup_tempfiles(temp_rds)
    }

    # Combine all
    counts_reprocessed <- bind_rows(counts_reprocessed_list)
    proportions_reprocessed <- bind_rows(proportions_reprocessed_list)
    tax_reprocessed <- bind_rows(tax_reprocessed_list)

    if (!file.exists(file.path(local, "rdp16classified.csv"))) {
      tax_reprocessed2 <- bind_rows(tax_reprocessed2_list)
      counts_reprocessed2 <- bind_rows(counts_reprocessed2_list)
      proportions_reprocessed2 <- bind_rows(proportions_reprocessed2_list)
      write.csv(tax_reprocessed2, file = file.path(local, "rdp16classified.csv"), row.names = TRUE)
    } else {
      counts_reprocessed2 <- bind_rows(counts_reprocessed2_list)
      proportions_reprocessed2 <- bind_rows(proportions_reprocessed2_list)
      tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
    }

  } else {
  # DELETE LATER WHEN REPROCESS HAS FINALIZED:
  counts_zip = file.path(local, "Marotz_2021_16S.csv.zip")
  counts_reprocessed = read_zipped_table(counts_zip) %>% as.data.frame()
  counts_reprocessed <- counts_reprocessed %>%
    rownames_to_column("SampleID") %>%
    mutate(SampleID = if_else(
      str_detect(SampleID, "raw"),
      SampleID,
      paste0(SampleID, ".PMA")
    ))
  counts_reprocessed <- counts_reprocessed %>% column_to_rownames("SampleID")
  if (!raw) {
    aligned = rename_and_align(counts_original = counts_reprocessed, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
    counts_reprocessed = aligned$counts_original
    original_names <- colnames(counts_reprocessed)
    counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
  }
  proportions_reprocessed <- sweep(counts_reprocessed, MARGIN = 1,STATS  = rowSums(counts_reprocessed), FUN = "/")
  tax_reprocessed <- data.frame(Taxa = colnames(counts_reprocessed), stringsAsFactors = FALSE)
  ############################################
  }

  if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_original = fill_na_zero_numeric(proportions_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
      proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
  }

  # ----- Return structured list -----
  return(list(
    counts = list(
      original = counts_original,
      reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
    ),
    tax = list(
      original = tax_original,
      reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
    ),
    scale = list(original = scale_metadata, reprocessed = scale_sra),
    metadata = list(original = metadata, reprocessed = sra)
  ))
}
