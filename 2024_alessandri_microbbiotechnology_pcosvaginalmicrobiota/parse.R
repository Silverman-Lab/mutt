parse_2024_alessandri_microbbiotechnology_pcosvaginalmicrobiota <- function(raw = FALSE, align = FALSE) {
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

  # ----- Local base directory -----
  local <- file.path("2024_alessandri_microbbiotechnology_pcosvaginalmicrobiota")

  # ----- File paths -----
  metadata_zip            <- file.path(local, "SraRunTable (1).csv.zip")
  motus_zip               <- file.path(local, "PRJNA1075117_motus_merged.tsv.zip")
  metaphlan4_zip          <- file.path(local, "PRJNA1075117_MetaPhlAn_merged.tsv.zip")
  scale_zip               <- file.path(local, "Alessandri2024_scale.csv.zip")
  metadataprocessed_zip   <- file.path(local, "Alessandri_2024_metadata.csv.zip")
  oldprocessed_zip        <- file.path(local, "Alessandri_2024_shotgunmetagenomics.csv.zip")

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

  # ----- Metadata -----
  sra          <- read_zipped_table(metadata_zip, row.names=NULL) %>% rename(Accession = Run, Sample = `Sample Name`)
  metadata <- read_zipped_table(metadataprocessed_zip, row.names=1) %>% rename(Accession = Run, Sample = `Sample.Name`)
  #no need to merge, just use metadata

  # ----- Scale -----
  scale <- read_zipped_table(scale_zip, row.names=1) %>% 
            mutate(log2_FC_Bacterial_load = ifelse(10^`Bacterial total count` > 0, log2(10^`Bacterial total count`),NA)) %>% 
            rename(log10_FC_Bacterial_load = `Bacterial total count`)
            

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
          aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, 
                            by_col="Sample", align = align, study_name=basename(local))
          df = aligned$reprocessed
          original_names <- colnames(df)
          df <- as.data.frame(lapply(df, as.numeric), row.names = rownames(df), col.names = original_names, check.names = FALSE)
        }
        cleanup_tempfiles(temp_dir)
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
    temp_dir <- tempfile("repro_")   
    dir.create(temp_dir)
    metaphlan4_files <- unzip(metaphlan4_zip, list = TRUE)
    metaphlan4_filename <- metaphlan4_files$Name[grepl("\\.tsv$", metaphlan4_files$Name)][1]
    if (!is.na(metaphlan4_filename)) {
        unzip(metaphlan4_zip, files = metaphlan4_filename, exdir = temp_dir, overwrite = TRUE)
        path <- file.path(temp_dir, metaphlan4_filename)
        df <- read_tsv(path)
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
        if (!raw) {
          aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, 
                            by_col="Sample", align = align, study_name=basename(local))
          df = aligned$reprocessed
          original_names <- colnames(df)
          df <- as.data.frame(lapply(df, as.numeric), row.names = rownames(df), col.names = original_names, check.names = FALSE)
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
    cleanup_tempfiles(temp_dir)
  }

  # ---- old processed, delete later ----
  counts_original <- read_zipped_table(oldprocessed_zip)
  if (!raw) {
    aligned = rename_and_align(counts_reprocessed = counts_original, metadata=metadata, scale=scale, 
                              by_col="Sample", align = align, study_name=basename(local))
    counts_original = aligned$reprocessed
    original_names <- colnames(counts_original)
    counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
  }
  proportions_original <- sweep(counts_original, MARGIN = 1,STATS  = rowSums(counts_original), FUN = "/")
  tax_original <- data.frame(Taxa = colnames(counts_original))

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
