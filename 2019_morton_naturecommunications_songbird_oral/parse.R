parse_2019_morton_naturecommunications_songbird_oral <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("stringr", "tidyverse", "matrixStats")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
            stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                ". Please install them before running this function.")
    }
    if (!is.logical(align)) {
      stop("align must be a logical value")
    }
    if (!is.logical(raw)) {
      stop("raw must be a logical value")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2019_morton_naturecommunications_songbird_oral")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJEB29169_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJEB29169_dada2_taxa.rds.zip")
    metadata_zip         <- file.path(local, "oral_trimmed_metadata.csv.zip")
    counts_zip           <- file.path(local, "2019_morton_songbird_oral_counts.RDS.zip")
    sra_zip              <- file.path(local, "SraRunTable (40).csv.zip")
    tax_zip              <- file.path(local, "taxonomy.tsv.zip")
    motus_zip            <- file.path(local, "PRJEB29169_motus_merged.tsv.zip")
    metaphlan4_zip       <- file.path(local, "PRJEB29169_MetaPhlAn_merged.tsv.zip")

    # ----- Initialize -----
    counts = NULL
    proportions = NULL
    tax = NULL
    metadata = NULL
    scale = NULL
    counts_reprocessed = NULL
    proportions_reprocessed = NULL
    tax_reprocessed = NULL
    tax_reprocessed2 = NULL
    counts_reprocessed2 = NULL
    proportions_reprocessed2 = NULL
    mOTU3_counts = NULL
    mOTU3_proportions = NULL
    mOTU3_tax = NULL
    MetaPhlAn4_counts = NULL
    MetaPhlAn4_proportions = NULL
    MetaPhlAn4_tax = NULL

    # ----- Metadata -----
    metadata <- read_zipped_table(metadata_zip, row.names=NULL)
    metadata <- metadata[, !is.na(names(metadata)) & names(metadata) != ""] 

    primarystudymetadata <- metadata %>%
      rename(Sample = SampleID) %>%
      mutate(across(c("flow cells/ul 1", "flow cells/ul 2"), as.numeric)) %>%
      mutate(across(c("qPCR cell/ul 1", "qPCR cell/ul 2", "qPCR cell/ul 3"), as.numeric)) %>%
      mutate(across(c("qPCR cell 5 min 1", "qPCR cell 5 min 2", "qPCR cell 5 min 3"), as.numeric)) %>%
      mutate(
        avg_FC_cells_ul = rowMeans(select(., "flow cells/ul 1", "flow cells/ul 2"), na.rm = TRUE),
        sd_FC_cells_ul = apply(select(., "flow cells/ul 1", "flow cells/ul 2"), 1, sd, na.rm = TRUE),
        avg_qpcr_cells_ul = rowMeans(select(., "qPCR cell/ul 1", "qPCR cell/ul 2", "qPCR cell/ul 3"), na.rm = TRUE),
        sd_qpcr_cells_ul = apply(select(., "qPCR cell/ul 1", "qPCR cell/ul 2", "qPCR cell/ul 3"), 1, sd, na.rm = TRUE),
        avg_qpcr_cells_5min = rowMeans(select(., "qPCR cell 5 min 1", "qPCR cell 5 min 2", "qPCR cell 5 min 3"), na.rm = TRUE),
        sd_qpcr_cells_5min = apply(select(., "qPCR cell 5 min 1", "qPCR cell 5 min 2", "qPCR cell 5 min 3"), 1, sd, na.rm = TRUE),
        SampleID = paste0(`participant-timepoint`, ".", Timepoint, ".", treatment, ".amplicon")
      ) 

    primarystudyscale = primarystudymetadata %>% select(SampleID, avg_FC_cells_ul, sd_FC_cells_ul, avg_qpcr_cells_5min, sd_qpcr_cells_5min, avg_qpcr_cells_ul, sd_qpcr_cells_ul) %>% 
      mutate(log2_FC_avg_cells_ul = ifelse(avg_FC_cells_ul > 0, log2(avg_FC_cells_ul), NA)) %>%
      mutate(log10_FC_avg_cells_ul = ifelse(avg_FC_cells_ul > 0, log10(avg_FC_cells_ul), NA)) %>%
      mutate(log2_FC_sd_cells_ul = ifelse(sd_FC_cells_ul > 0, log2(sd_FC_cells_ul), NA)) %>%
      mutate(log10_FC_sd_cells_ul = ifelse(sd_FC_cells_ul > 0, log10(sd_FC_cells_ul), NA)) %>%
      mutate(log2_qpcr_avg_cells_5min = ifelse(avg_qpcr_cells_5min > 0, log2(avg_qpcr_cells_5min), NA)) %>%
      mutate(log10_qpcr_avg_cells_5min = ifelse(avg_qpcr_cells_5min > 0, log10(avg_qpcr_cells_5min), NA)) %>%
      mutate(log2_qpcr_avg_cells_ul = ifelse(avg_qpcr_cells_ul > 0, log2(avg_qpcr_cells_ul), NA)) %>%
      mutate(log10_qpcr_avg_cells_ul = ifelse(avg_qpcr_cells_ul > 0, log10(avg_qpcr_cells_ul), NA)) %>%
      mutate(log2_qpcr_sd_cells_5min = ifelse(sd_qpcr_cells_5min > 0, log2(sd_qpcr_cells_5min), NA)) %>%
      mutate(log10_qpcr_sd_cells_5min = ifelse(sd_qpcr_cells_5min > 0, log10(sd_qpcr_cells_5min), NA)) %>%
      mutate(log2_qpcr_sd_cells_ul = ifelse(sd_qpcr_cells_ul > 0, log2(sd_qpcr_cells_ul), NA)) %>%
      mutate(log10_qpcr_sd_cells_ul = ifelse(sd_qpcr_cells_ul > 0, log10(sd_qpcr_cells_ul), NA))

    sra = read_zipped_table(sra_zip, row.names = NULL) %>% rename(Accession = Run)
    sra$SampleID <- paste0(sra$saliva_sample_id, ".", sra$timepoint, ".", sra$processing, ".", 
                          ifelse(sra$`Assay Type` == "WGS", "shotgun", "amplicon"))
    sra = sra %>%
      mutate(SampleID = case_when(
        grepl("^\\.\\d+\\..(shotgun|amplicon)$", SampleID) ~ NA_character_,
        grepl("^\\.NA\\..(shotgun|amplicon)$", SampleID) ~ NA_character_,
        TRUE ~ SampleID
      )) %>%
      mutate(SampleID = ifelse(is.na(SampleID), 
                              paste0(
                                sub("^([^.]+)\\.([^.]+)$", "\\1", `orig_name (exp)`), ".",
                                timepoint, ".",
                                sub("^([^.]+)\\.([^.]+)$", "\\2", `orig_name (exp)`), ".",
                                `plating (exp)`, ".",
                                ifelse(`Assay Type` == "WGS", "shotgun", "amplicon")
                              ),
                              SampleID)) %>%
      mutate(SampleID = case_when(
        grepl("^\\.\\d+\\...(shotgun|amplicon)$", SampleID) ~ NA_character_,
        grepl("^\\.NA\\...(shotgun|amplicon)$", SampleID) ~ NA_character_,
        TRUE ~ SampleID
      ))

    srascale <- sra %>% select(c("Accession", "fc_avg_cells_5_min", "fc_avg_cells_per_ul", "fc_cells_per_ul_r1", "fc_cells_per_ul_r2", "qpcr_median_16s_copies_per_2ul_dna", "all_flow_cells_5min_avg", "all_flow_cellsperul_avg", "all_qpcr_cells_5min_avg", "all_qpcr_cellsperul_avg", "live_flow_cells_5min_avg", "live_flow_cellsperul_avg", "live_qpcr_cells_5min_avg", "live_qpcr_cellsperul_avg")) %>%
      mutate(across(c("fc_avg_cells_5_min", "fc_avg_cells_per_ul", "fc_cells_per_ul_r1", "fc_cells_per_ul_r2", "qpcr_median_16s_copies_per_2ul_dna", "all_flow_cells_5min_avg", "all_flow_cellsperul_avg", "all_qpcr_cells_5min_avg", "all_qpcr_cellsperul_avg", "live_flow_cells_5min_avg", "live_flow_cellsperul_avg", "live_qpcr_cells_5min_avg", "live_qpcr_cellsperul_avg"), as.numeric)) %>%
      mutate(log2_fc_avg_cells_5_min = ifelse(fc_avg_cells_5_min > 0, log2(fc_avg_cells_5_min), NA)) %>%
      mutate(log10_fc_avg_cells_5_min = ifelse(fc_avg_cells_5_min > 0, log10(fc_avg_cells_5_min), NA)) %>%
      mutate(log2_fc_avg_cells_per_ul = ifelse(fc_avg_cells_per_ul > 0, log2(fc_avg_cells_per_ul), NA)) %>%
      mutate(log10_fc_avg_cells_per_ul = ifelse(fc_avg_cells_per_ul > 0, log10(fc_avg_cells_per_ul), NA)) %>%
      mutate(log2_fc_cells_per_ul_r1 = ifelse(fc_cells_per_ul_r1 > 0, log2(fc_cells_per_ul_r1), NA)) %>%
      mutate(log10_fc_cells_per_ul_r1 = ifelse(fc_cells_per_ul_r1 > 0, log10(fc_cells_per_ul_r1), NA)) %>%
      mutate(log2_fc_cells_per_ul_r2 = ifelse(fc_cells_per_ul_r2 > 0, log2(fc_cells_per_ul_r2), NA)) %>%
      mutate(log10_fc_cells_per_ul_r2 = ifelse(fc_cells_per_ul_r2 > 0, log10(fc_cells_per_ul_r2), NA)) %>%
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
      mutate(log10_live_flow_cellsperul_avg = ifelse(live_flow_cellsperul_avg > 0, log10(live_flow_cellsperul_avg), NA)) %>%
      mutate(log2_live_qpcr_cells_5min_avg = ifelse(live_qpcr_cells_5min_avg > 0, log2(live_qpcr_cells_5min_avg), NA)) %>%
      mutate(log10_live_qpcr_cells_5min_avg = ifelse(live_qpcr_cells_5min_avg > 0, log10(live_qpcr_cells_5min_avg), NA)) %>%
      mutate(log2_live_qpcr_cellsperul_avg = ifelse(live_qpcr_cellsperul_avg > 0, log2(live_qpcr_cellsperul_avg), NA)) %>%
      mutate(log10_live_qpcr_cellsperul_avg = ifelse(live_qpcr_cellsperul_avg > 0, log10(live_qpcr_cellsperul_avg), NA))

    ## Read Counts
    temp_rds <- tempfile("repro")
    dir.create(temp_rds)
    unzipped = unzip(counts_zip, exdir = dirname(temp_rds), overwrite = TRUE)
    counts_file <- unzipped[grep("2019_morton_songbird_oral_counts\\.RDS$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No 2019_morton_songbird_oral_counts.rds file found after unzip")
    counts <- as.data.frame(readRDS(counts_file)) %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Sample") 
    counts$Sample <- primarystudymetadata$SampleID[match(counts$Sample, primarystudymetadata$Sample)]
    counts <- counts %>% column_to_rownames(var = "Sample")
    ## Taxonomy Information
    raw_tax <- read_zipped_table(tax_zip, sep="\t", row.names = NULL)
    tax <- raw_tax %>%
        mutate(taxonomy = str_replace_all(Taxon, "\\s+", "")) %>%
        mutate(ogtaxonomy = taxonomy) %>%
        separate(
            taxonomy,
            into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
            sep = ";",
            extra = "drop",
            fill = "right"
        ) %>%
        mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                     ~gsub("^[kpcofgs]__", "", .))) %>%
        mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                     ~gsub("\\[|\\]|\\(|\\)", "", .))) %>% 
        mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                     ~gsub("_+", "_", .))) %>% 
        mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                     ~ifelse(. == "" | . == "__" | is.na(.), "unclassified", .)))
    tax = tax %>% rename(species = Species)
    tax = make_taxa_label(tax)
    tax = tax %>% rename(Species = species)
    tax <- tax[,c("Feature ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Taxa", "ogtaxonomy")]
    row.names(tax) <- tax$`Feature ID`

    if (!raw) {
        aligned = rename_and_align(counts_original = counts, metadata=primarystudymetadata, scale=primarystudyscale, by_col = "Sample", align = align, study_name = basename(local))
        counts <- aligned$counts_original
        matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
        colnames(counts) <- matched_taxa
        counts = collapse_duplicate_columns_exact(counts)
        original_names <- colnames(counts)
        counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
    }
    proportions <- sweep(counts, 1, rowSums(counts), FUN = "/")

    # ----- mOTU3 Reprocessed -----
    if (file.exists(motus_zip)) {
      # 1. create a private scratch folder
      temp_dir <- tempfile("motus_")
      dir.create(temp_dir)

      # 2. find the .tsv inside the ZIP
      motus_files    <- unzip(motus_zip, list = TRUE)
      motus_filename <- motus_files$Name[
                        grepl("\\.tsv$", motus_files$Name, ignore.case = TRUE)
                      ][1]

      if (!is.na(motus_filename)) {
        # 3. extract just that file and grab its full path
        unzipped    <- unzip(
                        motus_zip,
                        files     = motus_filename,
                        exdir     = temp_dir,
                        overwrite = TRUE
                      )
        motus_path  <- unzipped[1]

        # 4. read counts + set rownames
        df <- readr::read_tsv(motus_path, show_col_types = FALSE)
        rownames(df) <- df[[1]]
        df[[1]]      <- NULL

        # 5. optional alignment
        if (!raw) {
          aligned <- rename_and_align(
            counts_reprocessed = df,
            metadata          = srametadata,
            scale             = srascale,
            by_col            = "Accession",
            align             = align,
            study_name        = basename(local)
          )
          df <- aligned$reprocessed
        }

        # 6. numeric conversion + proportions
        df[]        <- lapply(df, as.numeric)
        proportions <- sweep(df, 1, rowSums(df), "/")

        # 7. simple taxonomy table from rownames
        tax_df <- tibble::tibble(taxa = rownames(df)) |>
          dplyr::mutate(taxa = stringr::str_trim(taxa)) |>
          tidyr::separate(
            taxa,
            into  = c("Kingdom","Phylum","Class","Order",
                      "Family","Genus","Species","Strain"),
            sep   = "\\s*;\\s*", extra = "drop", fill = "right"
          )
        rownames(tax_df) <- rownames(df)

        # 8. assign out
        mOTU3_counts       <- df
        mOTU3_proportions  <- proportions
        mOTU3_tax          <- tax_df
      }

      # 9. clean up only your private folder
      cleanup_tempfiles(temp_dir)
    }


    # ----- MetaPhlAn4 Reprocessed -----
    if (file.exists(metaphlan4_zip)) {
      # 1. private scratch folder
      temp_dir <- tempfile("mp4_")
      dir.create(temp_dir)

      # 2. locate the .tsv in the ZIP
      mp4_files    <- unzip(metaphlan4_zip, list = TRUE)
      mp4_filename <- mp4_files$Name[
                        grepl("\\.tsv$", mp4_files$Name, ignore.case = TRUE)
                      ][1]

      if (!is.na(mp4_filename)) {
        # 3. extract and capture full path
        unzipped  <- unzip(
                      metaphlan4_zip,
                      files     = mp4_filename,
                      exdir     = temp_dir,
                      overwrite = TRUE
                    )
        path      <- unzipped[1]

        # 4. read + set rownames
        df <- readr::read_tsv(path, show_col_types = FALSE)
        rownames(df) <- df[[1]]
        df[[1]]      <- NULL

        # 5. optional alignment
        if (!raw) {
          aligned <- rename_and_align(
            counts_reprocessed = df,
            metadata          = srametadata,
            scale             = srascale,
            by_col            = "Accession",
            align             = align,
            study_name        = basename(local)
          )
          df <- aligned$reprocessed
        }

        # 6. numeric + proportions
        df[]        <- lapply(df, as.numeric)
        proportions <- sweep(df, 1, rowSums(df), "/")

        # 7. taxonomy table
        tax_df <- tibble::tibble(taxa = rownames(df)) |>
          dplyr::mutate(taxa = stringr::str_trim(taxa)) |>
          tidyr::separate(
            taxa,
            into  = c("Kingdom","Phylum","Class","Order",
                      "Family","Genus","Species","Strain"),
            sep   = "\\s*;\\s*", extra = "drop", fill = "right"
          )
        rownames(tax_df) <- rownames(df)

        # 8. assign out
        MetaPhlAn4_counts      <- df
        MetaPhlAn4_proportions <- proportions
        MetaPhlAn4_tax         <- tax_df
      }

      # 9. tidy up
      cleanup_tempfiles(temp_dir)
    }
    
    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
        unzipped = unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)
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
        unzipped = unzip(repro_tax_zip, exdir = dirname(temp_rds), overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))
        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=srametadata, scale=srascale, by_col = "Accession", align = align, study_name = basename(local))
            counts_reprocessed <- aligned$reprocessed
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
    }

    if (!raw) {
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        counts = fill_na_zero_numeric(counts)
        proportions = fill_na_zero_numeric(proportions)
        counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
        proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    }   

    cleanup_tempfiles(temp_rds)

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts,
            reprocessed = list(amplicon = list(
              rdp19 = counts_reprocessed, 
              rdp16 = counts_reprocessed2
            ),
            shotgun = list(
              mOTU3 = mOTU3_counts,
              MetaPhlAn4 = MetaPhlAn4_counts
            )
            )
        ),
        proportions = list(
            original = proportions,
            reprocessed = list(amplicon = list(
              rdp19 = proportions_reprocessed, 
              rdp16 = proportions_reprocessed2
            ),
            shotgun = list(
              mOTU3 = mOTU3_proportions,
              MetaPhlAn4 = MetaPhlAn4_proportions
            )
            )
        ),
        tax = list(
            original = tax,
            reprocessed = list(amplicon = list(
              rdp19 = tax_reprocessed, 
              rdp16 = tax_reprocessed2
            ),
            shotgun = list(
              mOTU3 = mOTU3_tax,
              MetaPhlAn4 = MetaPhlAn4_tax
            )
            )
        ),
        scale = list(original = primarystudyscale, 
                    reprocessed = srascale
                    ),
        metadata = list(original = primarystudymetadata, 
                    reprocessed = sra
                    )
    ))
}

