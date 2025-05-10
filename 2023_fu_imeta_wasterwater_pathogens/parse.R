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
    #     proportions <- as.data.frame(t(rowsum(t(proportions), group = colnames(proportions))))
    # }

  mockproportions = read_zipped_table(mock_prop_zip, row.names=NULL) %>% as.data.frame()
  mocktax <- tibble(taxonomy = colnames(mockproportions))

  # ----- Initialize everything as NA -----
  mOTU3_counts <- NA
  mOTU3_proportions <- NA
  mOTU3_tax <- NA

  MetaPhlAn4_counts <- NA
  MetaPhlAn4_proportions <- NA
  MetaPhlAn4_tax <- NA

  # ----- mOTU3 Reprocessed -----
  if (file.exists(motus_zip)) {
    motus_files <- unzip(motus_zip, list = TRUE)
    motus_filename <- motus_files$Name[grepl("\\.tsv$", motus_files$Name)][1]
    if (!is.na(motus_filename)) {
        temp_dir <- tempdir()
        unzip(motus_zip, files = motus_filename, exdir = temp_dir, overwrite = TRUE)
        motus_path <- file.path(temp_dir, motus_filename)
        df <- read_tsv(motus_path)
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
            df = aligned$reprocessed
        }
        proportions <- sweep(df, 1, rowSums(df), FUN = "/")
        tax_df <- data.frame(taxa = rownames(df)) %>%
        mutate(taxa = str_trim(taxa)) %>%
        separate(taxa,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                sep = "\\s*;\\s*", extra = "drop", fill = "right")
        rownames(tax_df) <- rownames(df)

        mOTU3_counts <- df
        mOTU3_proportions <- proportions
        mOTU3_tax <- tax_df
    }
  }

  # ----- MetaPhlAn4 Reprocessed -----
  if (file.exists(metaphlan4_zip)) {
    metaphlan4_files <- unzip(metaphlan4_zip, list = TRUE)
    metaphlan4_filename <- metaphlan4_files$Name[grepl("\\.tsv$", metaphlan4_files$Name)][1]
    if (!is.na(metaphlan4_filename)) {
        temp_dir <- tempdir()
        unzip(metaphlan4_zip, files = metaphlan4_filename, exdir = temp_dir, overwrite = TRUE)
        path <- file.path(temp_dir, metaphlan4_filename)
        df <- read_tsv(path)
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
            df = aligned$reprocessed
        }
        proportions <- sweep(df, 1, rowSums(df), FUN = "/")
        tax_df <- data.frame(taxa = rownames(df)) %>%
        mutate(taxa = str_trim(taxa)) %>%
        separate(taxa,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                sep = "\\s*;\\s*", extra = "drop", fill = "right")
        rownames(tax_df) <- rownames(df)

        MetaPhlAn4_counts <- df
        MetaPhlAn4_proportions <- proportions
        MetaPhlAn4_tax <- tax_df
    }
  }


  if (!raw) {
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
      proportions = fill_na_zero_numeric(proportions)
      mockproportions = fill_na_zero_numeric(mockproportions)
  }

  return(list(
    counts      = list(original = list(
                        mocks = NA, 
                        samples = NA
                        ),
                      reprocessed = list(
                        mock = NA,
                        samples = counts_reprocessed
                        )
                      ),
    proportions = list(original = list(
                        mock = mockproportions, 
                        samples = proportions # NEEDS TO BE FIXED because bins do not match sra
                      ),
                      reprocessed = list(
                        mock = NA,
                        samples = proportions_reprocessed 
                      )
                      ),
    tax         = list(original = list(
                        mock = mocktax, 
                        samples = tax
                        ),
                      reprocessed = list(
                        mock = NA,
                        samples = tax_reprocessed
                        ),
                      ),
    scale       = list(
                    mock = scale,
                    samples = NA #IDK
                      ), 
    metadata    = list(
                    mock = mockmetadata,
                    samples = sra
                      )
  ))
}