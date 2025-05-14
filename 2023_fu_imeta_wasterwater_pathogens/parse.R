parse_2023_fu_imeta_wasterwater_pathogens <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "R.utils")
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
  local <- file.path("2023_fu_imeta_wasterwater_pathogens")

  # ----- File paths -----
  mock_scale_zip       <- file.path(local, "mock_scale.csv.zip")
  mock_prop_zip        <- file.path(local, "sheet1_pivoted.csv.zip")
  samples_prop_zip     <- file.path(local, "sheet_s5_cleaned.csv.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  motus_zip            <- file.path(local, "PRJNA860773_motus_merged.tsv.zip")
  metaphlan4_zip       <- file.path(local, "PRJNA860773_MetaPhlAn_merged.tsv.zip")

  counts = NA
  proportions = NA
  tax = NA
  scale = NA
  metadata = NA
  sra = NA
  mock = NA

  # scale
  mockscale = read_zipped_table(mock_scale_zip, row.names=NULL) %>% as.data.frame() %>% 
                      mutate(log2_FC_mean = ifelse(Flow_mean > 0, log2(Flow_mean), NA)) %>%
                      mutate(log2_FC_SD = ifelse(Flow_SD > 0, log2(Flow_SD), NA)) %>%
                      mutate(log2_FC_qPCR_mean = ifelse(qPCR_mean > 0, log2(qPCR_mean), NA)) %>%
                      mutate(log2_FC_qPCR_SD = ifelse(qPCR_SD > 0, log2(qPCR_SD), NA)) %>% 
                      mutate(log10_FC_mean = ifelse(Flow_mean > 0, log10(Flow_mean), NA)) %>%
                      mutate(log10_FC_SD = ifelse(Flow_SD > 0, log10(Flow_SD), NA)) %>%
                      mutate(log10_FC_qPCR_mean = ifelse(qPCR_mean > 0, log10(qPCR_mean), NA)) %>%
                      mutate(log10_FC_qPCR_SD = ifelse(qPCR_SD > 0, log10(qPCR_SD), NA))


  # metadata
  sra = read_zipped_table(metadata_zip, row.names=NULL) %>% as.data.frame() %>% rename(Accession = Run)
  # proportions

  proportions = read_zipped_table(samples_prop_zip, row.names=NULL) %>% as.data.frame()
  proportions <- proportions %>%
  mutate(`Relative abundance %` = as.numeric(`Relative abundance %`)) %>%
  filter(!is.na(`Relative abundance %`), !is.na(classification), classification != "")  

  # Pivot into wide format: rows are bins, columns are classifications
  proportions <- proportions %>%
    select(bins, classification, `Relative abundance %`) %>%
    pivot_wider(
      names_from = classification,
      values_from = `Relative abundance %`,
      values_fill = 0
    ) %>% 
    column_to_rownames("bins")
  


  # tax 
  tax <- tibble(
    taxonomy = colnames(proportions)
  )

  tax$ogtaxonomy <- tax$taxonomy
    tax <- tax %>%
    separate(taxonomy,
            into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
            sep = ";", fill = "right") %>%
    mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))
    tax <- tax %>%
        mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), NA, .)))
    tax[is.na(tax)] <- "unclassified"
    tax = make_taxa_label(tax)
    rownames(tax) <- tax$ogtaxonomy
    # if (!raw) {
    #     matched_taxa <- tax$Taxa[match(colnames(proportions), rownames(tax))]
    #     colnames(proportions) <- matched_taxa
    #     proportions = collapse_duplicate_columns_exact(proportions)
    #     original_names <- colnames(proportions)
    #     proportions <- as.data.frame(lapply(proportions, as.numeric), row.names = rownames(proportions), col.names = original_names, check.names = FALSE)
    # }

  mockproportions = read_zipped_table(mock_prop_zip, row.names=1) %>% as.data.frame()
  mocktax <- tibble(taxonomy = colnames(mockproportions))
  if (!raw) {
    aligned <- rename_and_align(proportions_original = mockproportions, metadata = metadata, scale = mockscale, by_col = "Sample", align = align, study_name = basename(local))
    mockproportions = aligned$proportions_original
    original_names <- colnames(mockproportions)
    mockproportions <- as.data.frame(lapply(mockproportions, as.numeric), row.names = rownames(mockproportions), col.names = original_names, check.names = FALSE)
  }


  # ----- Initialize everything as NA -----
  mOTU3_counts <- NA
  mOTU3_proportions <- NA
  mOTU3_tax <- NA

  MetaPhlAn4_counts <- NA
  MetaPhlAn4_proportions <- NA
  MetaPhlAn4_tax <- NA
  
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
          by_col            = "Sample",
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
          metadata          = metadata,
          scale             = scale,
          by_col            = "Sample",
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


  if (!raw) {

      mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
      mOTU3_proportions = fill_na_zero_numeric(mOTU3_proportions)
      MetaPhlAn4_counts = fill_na_zero_numeric(MetaPhlAn4_counts)
      MetaPhlAn4_proportions = fill_na_zero_numeric(MetaPhlAn4_proportions)
      mockproportions = fill_na_zero_numeric(mockproportions)
      proportions = fill_na_zero_numeric(proportions)
  }

  return(list(
    counts      = list(original = list(
                        mocks = NA, 
                        samples = NA
                        ),
                      reprocessed = list(
                        mock = NA,
                        samples = list(mOTU3 = mOTU3_counts, MetaPhlAn4 = MetaPhlAn4_counts)
                        )
                      ),
    proportions = list(original = list(
                        mock = mockproportions, 
                        samples = proportions # NEEDS TO BE FIXED because bins do not match sra
                      ),
                      reprocessed = list(
                        mock = NA,
                        samples = list(mOTU3 = mOTU3_proportions, MetaPhlAn4 = MetaPhlAn4_proportions)
                      )
                      ),
    tax         = list(original = list(
                        mock = mocktax, 
                        samples = tax
                        ),
                      reprocessed = list(
                        mock = NA,
                        samples = list(mOTU3 = mOTU3_tax, MetaPhlAn4 = MetaPhlAn4_tax)
                        )
                      ),
    scale       = list(
                    mock = mockscale,
                    samples = scale #IDK
                      ), 
    metadata    = list(
                    mock = NA,
                    samples = sra
                      )
  ))
}