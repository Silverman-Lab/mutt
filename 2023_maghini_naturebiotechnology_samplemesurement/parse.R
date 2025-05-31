parse_2023_maghini_naturebiotechnology_samplemesurement <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function. For taxizedb, you can use:
           install.packages('remotes')
           remotes::install_github('ropensci/taxizedb')
         ")
  }
  library(tibble)
  library(tidyverse)
  library(readr)
  #library(taxizedb)

  # ----- Local base directory -----
  local <- file.path("2023_maghini_naturebiotechnology_samplemesurement")

  # ----- File paths -----
  motus_zip      <- file.path(local, "PRJNA940499_motus_merged.tsv.zip")
  metaphlan4_zip <- file.path(local, "PRJNA940499_MetaPhlAn_merged_counts.tsv.zip")
  scale_zip      <- file.path(local, "Maghini2023_scale.csv.zip")
  metadata_zip   <- file.path(local, "Maghini_2023_metadata.csv.zip")
  original_zip   <- file.path(local, "Maghini_2023_shotgunmetagenomics.csv.zip")
  metadatafromsra_path <- file.path(local, "SraRunTable (3).csv.zip")

  # ----- Initialize everything as NA -----
  counts_original <- NA
  proportions_original <- NA
  tax_original <- NA
  mOTU3_counts <- NA
  mOTU3_proportions <- NA
  mOTU3_tax <- NA
  MetaPhlAn4_counts <- NA
  MetaPhlAn4_proportions <- NA
  MetaPhlAn4_tax <- NA

  # ----- Read and filter qPCR data -----
  plate_list <- list.files(path = local, pattern = "qPCR_plate", full.names = TRUE)
  qPCR <- plate_list %>%
    map_dfr(read_csv) %>%
    mutate(ID = gsub("_", "-", SampleName)) %>%
    dplyr::select(Plate, ID, Donor, Condition, PCR_Replicate, Replicate, logCopyNumber, CopyNumber) 

  scale <- qPCR %>%
    select(ID, logCopyNumber, CopyNumber) %>%
    group_by(ID) %>%
    summarise(
      mean_CopyNumber    = mean(CopyNumber, na.rm = TRUE),
      sd_CopyNumber      = sd(CopyNumber, na.rm = TRUE),
      mean_logCopyNumber = mean(logCopyNumber, na.rm = TRUE),
      sd_logCopyNumber   = sd(logCopyNumber, na.rm = TRUE),
      n_reps             = n(),
      .groups = "drop"
    ) %>%
  mutate(
    log2_CopyNumber = ifelse(mean_CopyNumber > 0, log2(mean_CopyNumber), NA),
    log10_CopyNumber = ifelse(mean_CopyNumber > 0, log10(mean_CopyNumber), NA)
  )

  # ----- Extract and deduplicate metadata -----
  metadata <- qPCR %>%
    dplyr::select(ID, Donor, Condition, Replicate) %>%
    distinct()

  # ----- SRA metadata and merge -----
  metadatafromsra <- read_zipped_table(metadatafromsra_path, row.names=NULL) %>%
    mutate(ID = gsub("(_DNA)", "", `Library Name`)) %>% rename(Accession = Run)

  scale$ID <- gsub("-", "_", scale$ID)
  metadata$ID <- gsub("-", "_", metadata$ID)

  metadata <- full_join(
    metadatafromsra, metadata,
    by = "ID"
  )

  metadata <- metadata %>% mutate(
    across(c(`Assay Type`, BioProject, BioSampleModel, Consent,
             `DATASTORE filetype`, `DATASTORE provider`, `DATASTORE region`,
             geo_loc_name_country, geo_loc_name_country_continent,
             Instrument, isolation_source, LibraryLayout, LibrarySelection,
             LibrarySource, Organism, Platform),
           ~ factor(na_if(., ""))),
    
    Donor = factor(Donor.x),
    Condition = factor(Condition),
    Replicate = factor(Replicate),
    ReleaseDate = as.POSIXct(ReleaseDate, format = "%Y-%m-%dT%H:%M:%SZ"),
    create_date = as.POSIXct(create_date, format = "%Y-%m-%dT%H:%M:%SZ"),
    Collection_Date = as.Date(paste0(Collection_Date, "-01-01")),
    across(c(Accession, BioSample, Experiment, `Sample Name`), factor),
    across(c(AvgSpotLen, Bases, Bytes), as.numeric),
    samp_collect_device = factor(replace_na(na_if(samp_collect_device, ""), "No preservative")),
    samp_mat_process = ordered(case_when(
      grepl("-80", samp_mat_process) ~ "-80C",
      grepl("23", samp_mat_process) ~ "23C",
      grepl("40", samp_mat_process) ~ "40C",
      TRUE ~ NA_character_
    ), levels = c("-80C", "23C", "40C"))
  ) %>% select(-Donor.y, -Donor.x) 

  ## vvvvvvvvv UNCOMMENT WHEN FIXED vvvvvvvvv


  # ----- Read taxonomic counts -----

  # tax_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  # read_bracken_table <- function(tax_level) {
  #   zip_path <- file.path(local, paste0("bracken_", tax_level, "_reads.txt.zip"))
  #   txt_path <- file.path(local, paste0("bracken_", tax_level, "_reads.txt"))
  #   if (file.exists(zip_path)) {
  #     tmp_dir <- tempdir()
  #     extracted <- tryCatch(unzip(zip_path, exdir = tmp_dir), error = function(e) NULL)
  #     txt_file <- extracted[grepl("\\.txt$", extracted)][1]
  #     if (!is.null(txt_file) && file.exists(txt_file)) {
  #       return(read.table(txt_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>% t() %>% as.data.frame())
  #     } else {
  #       warning("No .txt file found inside: ", zip_path)
  #       return(NULL)
  #     }
  #   }
  #   if (file.exists(txt_path)) {
  #     return(read.table(txt_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>% t() %>% as.data.frame())
  #   }

  #   warning("Neither ZIP nor TXT file found for level: ", tax_level)
  #   return(NULL)
  # }
  # raw_counts <- setNames(lapply(tax_levels, read_bracken_table), tax_levels)

  # renamed <- imap(raw_counts, function(df, level) {
  #   orig <- colnames(df)
  #   list(
  #     counts = df,
  #     tax    = tibble(taxonomy = orig)
  #   )
  # })
         
  # counts_original  <- map(renamed, "counts")
  # tax_original     <- map(renamed, "tax")

  # counts_original <- counts_original$species

  # species_names <- tax_original$species$taxonomy
  # taxizedb::db_download_ncbi()
  # file.copy(taxizedb::db_path("ncbi"), "ncbi_taxonomy.sqlite")
  # lineages <- classification(species_names, db = "ncbi", dbpath = "ncbi_taxonomy.sqlite")
  # standard_ranks <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  # tax <- imap_dfr(lineages, function(df, taxon) {
  #   if (is.null(df)) return(NULL)
  #   df %>%
  #     filter(rank %in% standard_ranks) %>%
  #     distinct(rank, .keep_all = TRUE) %>%
  #     pivot_wider(names_from = rank, values_from = name) %>%
  #     mutate(Species = taxon)
  # }, .id = "lookup_name") %>%
  #   rename(Kingdom = superkingdom)
  # tax = make_taxa_label(tax)
  # rownames(tax) <- tax$Species
  # if (!raw) {
  #     aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="ID", align = align, study_name=basename(local))
  #     counts_original = aligned$counts_original
  #     matched_taxa <- tax$Taxa[match(colnames(counts_original), rownames(tax))]
  #     colnames(counts_original) <- matched_taxa
  #     counts_original <- collapse_duplicate_columns_exact(counts_original)
  #     original_names <- colnames(counts_original)
  #     counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
  # }

  # # ------- Proportions ----------
  # proportions_original <- map(counts_original, function(df) {
  #   prop  <- sweep(df, 1, rowSums(df), FUN = "/")
  #   return(as_tibble(prop))
  # })

  # cleanup_tempfiles(temp_dir)

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
          metadata          = metadata,
          scale             = scale,
          by_col            = "ID",
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
      df <- readr::read_tsv(path, show_col_types = FALSE) %>% as.data.frame() %>% column_to_rownames("clade") %>% t() %>% as.data.frame()

      # 5. optional alignment
      if (!raw) {
      aligned <- rename_and_align(
          counts_reprocessed = df,
          metadata          = metadata,
          scale             = scale,
          by_col            = "ID",
          align             = align,
          study_name        = basename(local)
      )
      df <- aligned$reprocessed
      }

      # 6. numeric + proportions
      df[]        <- lapply(df, as.numeric)
      proportions <- sweep(df, 1, rowSums(df), "/")

      # 7. taxonomy table
      tax_df <- data.frame(taxa = colnames(df)) %>%
      mutate(taxa = str_trim(taxa)) %>%
      separate(
          taxa,
          into  = c("Kingdom","Phylum","Class","Order",
                  "Family","Genus","Species","Strain"),
          sep   = "\\|",
          extra = "drop",
          fill  = "right"
      ) %>%
      # remove the leading letter__ (e.g. "k__", "p__") from every column
      mutate(across(Kingdom:Strain, ~ str_remove(.x, "^[a-z]__")))
      rownames(tax_df) <- colnames(df)
      tax_df = make_taxa_label(tax_df)

      # 8. assign out
      MetaPhlAn4_counts      <- df
      MetaPhlAn4_proportions <- proportions
      MetaPhlAn4_tax         <- tax_df
  }

  # 9. tidy up
  cleanup_tempfiles(temp_dir)
  }

  # Can delete later:
  counts_original = read_zipped_table(original_zip)
  tax_original <- data.frame(Taxa = colnames(counts_original))
  if (!raw) {
    aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="ID", align = align, study_name=basename(local))
    counts_original = aligned$counts_original
    original_names <- colnames(counts_original)
    counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
  }
  proportions_original <- sweep(counts_original, MARGIN = 1, STATS = rowSums(counts_original), FUN = "/")
  
  
  # Dont delete:
  if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
      proportions_original = fill_na_zero_numeric(proportions_original)
      MetaPhlAn4_counts = fill_na_zero_numeric(MetaPhlAn4_counts)
      mOTU3_proportions = fill_na_zero_numeric(mOTU3_proportions)
      MetaPhlAn4_proportions = fill_na_zero_numeric(MetaPhlAn4_proportions)
  }

  # ----- Return -----
  return(list(
    counts = list(
      original = counts_original,
      reprocessed = list(
        mOTU3 = mOTU3_counts,
        MetaPhlAn4 = MetaPhlAn4_counts
      )
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = list(
        mOTU3 = mOTU3_proportions,
        MetaPhlAn4 = MetaPhlAn4_proportions
      )
    ),
    tax = list(
      original = tax_original,
      reprocessed = list(
        mOTU3 = mOTU3_tax,
        MetaPhlAn4 = MetaPhlAn4_tax
      )
    ),
    scale = scale,
    metadata = metadata
  ))
}
