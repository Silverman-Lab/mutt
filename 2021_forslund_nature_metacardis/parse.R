parse_2021_forslund_nature_metacardis <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
      stop(
          "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
          ". Please install them before running this function."
      )
  }
  if (!is.logical(raw) || length(raw) != 1) {
    stop("`raw` must be a single logical value (TRUE or FALSE)")
  }
  if (!is.logical(align) || length(align) != 1) {
    stop("`align` must be a single logical value (TRUE or FALSE)")
  }

  # Load needed libraries
  library(tidyverse)
  library(readxl)
  library(readr)

  # -------- local path -----------
  local                       <- file.path("2021_forslund_nature_metacardis")

  # -------- file paths -----------
  counts_zip                  <- file.path(local, "MetaCardis_mOTUs_v25.tsv.zip")
  metamatmetformin_zip        <- file.path(local, "metaMatMetformin.RData.zip")
  scale_file                  <- file.path(local, "MetaCardis_load.tsv.zip")
  metadata_zip                <- file.path(local, "41586_2021_4177_MOESM3_ESM.xlsx.zip")
  sra_zips                    <- c(
    file.path(local, "SraRunTable (36).csv.zip"),
    file.path(local, "SraRunTable (35).csv.zip"),
    file.path(local, "SraRunTable (34).csv.zip")
  )
  repro_motus_zips            <- c(
    file.path(local, "PRJEB41311_motus_merged.tsv.zip"), # ALSO INCLUDED IN FROMENTIN 2022
    file.path(local, "PRJEB38742_motus_merged.tsv.zip"), # ALSO INCLUDED IN FROMENTIN 2022
    file.path(local, "PRJEB37249_motus_merged.tsv.zip")  # ALSO INCLUDED IN 2020_vierasilva_nature_BMIS
  )

  repro_metaphlan_zips        <- c(
    file.path(local, "PRJEB41311_MetaPhlAn_merged.tsv.zip"), # ALSO INCLUDED IN FROMENTIN 2022
    file.path(local, "PRJEB38742_MetaPhlAn_merged.tsv.zip"), # ALSO INCLUDED IN FROMENTIN 2022
    file.path(local, "PRJEB37249_MetaPhlAn_merged.tsv.zip")  # ALSO INCLUDED IN 2020_vierasilva_nature_BMIS
  )
  

  # reduced_feature_zip <- paste0(local, "reduced_feature.RData.zip")
  # if (file.exists(reduced_feature_zip)) {
  #   message("Extracting and loading reduced_feature data from ", reduced_feature_zip)
    
  #   temp_dir <- tempdir()
  #   unzip(reduced_feature_zip, exdir = temp_dir)
  #   rdata_files <- list.files(temp_dir, pattern = "\\.RData$", full.names = TRUE)
    
  #   for (rdata_file in rdata_files) {
  #     env <- new.env()
  #     load(rdata_file, envir = env)

  #     if (exists("reduced_feature", envir = env)) {
  #       out$reduced_feature <- env$reduced_feature
  #     } else {
  #       warning("reduced_feature not found in ", rdata_file)
  #     }
  #   }
  # } else {
  #   warning("Reduced_feature zip file not found: ", reduced_feature_zip)
  # }
  
  # ------- metadata -----------

  # if (file.exists(metamatmetformin_zip)) {
  #   message("Extracting and loading metamatmetformin data from ", metamatmetformin_zip)
    
  #   temp_dir <- tempdir()
  #   unzip(metamatmetformin_zip, exdir = temp_dir)
  #   rdata_files <- list.files(temp_dir, pattern = "\\.RData$", full.names = TRUE)
    
  #   for (rdata_file in rdata_files) {
  #     env <- new.env()
  #     load(rdata_file, envir = env)

  #     if (exists("metamatmetformin", envir = env)) {
  #       metadata_small <- env$metamatmetformin
  #     } else {
  #       warning("metamatmetformin not found in ", rdata_file)
  #     }
  #   }
  # } else {
  #   warning("Metamatmetformin zip file not found: ", metamatmetformin_zip)
  # }

  for (sra_zip in sra_zips) {
    if (file.exists(sra_zip)) {
      zinfo <- unzip(sra_zip, list = TRUE)
      csv_name <- zinfo$Name[grepl("\\.csv$", zinfo$Name)][1]
      if (!is.na(csv_name)) {
        tmpdir <- tempfile("sra_csv_")
        dir.create(tmpdir)
        unzip(sra_zip, files = csv_name, exdir = tmpdir)
        csv_path <- file.path(tmpdir, csv_name)

        sra_df <- read_csv(csv_path, show_col_types = FALSE)
        sra_list[[length(sra_list) + 1]] <- sra_df

        unlink(tmpdir, recursive = TRUE)
      }
    }
  }

  combined_sra_df <- bind_rows(sra_list)

  if (file.exists(metadata_zip)) {
    temp_dir <- tempdir()
    unzip(metadata_zip, exdir = temp_dir)
    excel_files <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)

    if (length(excel_files) > 0) {
      metadata_file <- excel_files[1]

      metadata_df <- read_excel(metadata_file, sheet = 3) %>%
        mutate(
          SampleID = if_else(
            grepl("^x", SampleID, ignore.case = FALSE),
            paste0("M0", SampleID),
            SampleID
          )
        )
    }

    unlink(tmpdir, recursive = TRUE)
  } else {
      warning("No .xlsx file found in ", metadata_zip)
  }

  metadata <- full_join(
    combined_sra_df,
    metadata_df,
    by = c("Submitter_Id" = "SampleID") 
  )

  # ----- scale --------
  if (file.exists(scale_file)) {
    scale_df <- read.table(unz(scale_file, "MetaCardis_load.tsv"), 
                            header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    scale <- as.matrix(scale_df)
  } else {
  warning("Scale file not found: ", scale_file)
  }

    # ----- original counts ---------
  if (file.exists(counts_zip)) {
    counts_df <- read.table(unz(counts_zip, "MetaCardis_mOTUs_v25.tsv"), 
                            header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    counts <- as.matrix(counts_df)
    if (!raw) {
      align <- rename_and_align(counts_original = counts, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
      counts <- align$counts_original
    }
    row_sums <- rowSums(counts)
    proportions <- sweep(counts, 1, row_sums, FUN = "/")
    proportions[is.nan(proportions)] <- 0  

  } else {
  warning("Counts file not found: ", counts_zip)
  }

  # Initialize empty dataframes
  mOTU3_counts <- NULL
  mOTU3_proportions <- NULL
  mOTU3_tax <- NULL
  MetaPhlAn4_counts <- NULL
  MetaPhlAn4_proportions <- NULL
  MetaPhlAn4_tax <- NULL

  # ----- mOTU3 Reprocessed -----
  for (motus_zip in repro_motus_zips) {
    if (file.exists(motus_zip)) {
      motus_files <- unzip(motus_zip, list = TRUE)
      motus_filename <- motus_files$Name[grepl("\\.tsv$", motus_files$Name)][1]
      if (!is.na(motus_filename)) {
        temp_dir <- tempdir()
        unzip(motus_zip, files = motus_filename, exdir = temp_dir, overwrite = TRUE)
        motus_path <- file.path(temp_dir, motus_filename)
        df <- read_tsv(motus_path, show_col_types = FALSE)
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
        if (!raw) {
          align <- rename_and_align(counts_reprocessed = df, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
          df <- align$reprocessed
        }

        # Normalize to proportions
        prop <- apply(df, 2, function(col) col / sum(col))

        # Taxonomy table
        tax_df <- data.frame(taxa = rownames(df)) %>%
          mutate(taxa = str_trim(taxa)) %>%
          separate(taxa,
                  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                  sep = "\\s*;\\s*", extra = "drop", fill = "right")
        rownames(tax_df) <- rownames(df)

        # Merge across studies
        if (is.null(mOTU3_counts)) {
          mOTU3_counts      <- df
          mOTU3_proportions <- prop
          mOTU3_tax         <- tax_df
        } else {
          mOTU3_counts      <- full_join(mOTU3_counts, df, by = "row.names") %>%
                              column_to_rownames(var = "Row.names")
          mOTU3_proportions <- full_join(mOTU3_proportions, prop, by = "row.names") %>%
                              column_to_rownames(var = "Row.names")
          mOTU3_tax         <- full_join(mOTU3_tax, tax_df, by = "row.names") %>%
                              column_to_rownames(var = "Row.names")
        }
      }
    }
  }

  # ----- MetaPhlAn4 Reprocessed -----
  for (metaphlan4_zip in repro_metaphlan_zips) {
    if (file.exists(metaphlan4_zip)) {
      metaphlan4_files <- unzip(metaphlan4_zip, list = TRUE)
      metaphlan4_filename <- metaphlan4_files$Name[grepl("\\.tsv$", metaphlan4_files$Name)][1]
      if (!is.na(metaphlan4_filename)) {
        temp_dir <- tempdir()
        unzip(metaphlan4_zip, files = metaphlan4_filename, exdir = temp_dir, overwrite = TRUE)
        path <- file.path(temp_dir, metaphlan4_filename)
        df <- read_tsv(path, show_col_types = FALSE)
        rownames(df) <- df[[1]]
        df[[1]] <- NULL
        if (!raw) {
          align <- rename_and_align(counts_reprocessed = df, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
          df <- align$reprocessed
        }

        # Normalize to proportions
        prop <- apply(df, 2, function(col) col / sum(col))

        # Taxonomy table
        tax_df <- data.frame(taxa = rownames(df)) %>%
          mutate(taxa = str_trim(taxa)) %>%
          separate(taxa,
                  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
                  sep = "\\s*;\\s*", extra = "drop", fill = "right")
        rownames(tax_df) <- rownames(df)

        # Merge across studies
        if (is.null(MetaPhlAn4_counts)) {
          MetaPhlAn4_counts      <- df
          MetaPhlAn4_proportions <- prop
          MetaPhlAn4_tax         <- tax_df
        } else {
          MetaPhlAn4_counts      <- full_join(MetaPhlAn4_counts, df, by = "row.names") %>%
                                    column_to_rownames(var = "Row.names")
          MetaPhlAn4_proportions <- full_join(MetaPhlAn4_proportions, prop, by = "row.names") %>%
                                    column_to_rownames(var = "Row.names")
          MetaPhlAn4_tax         <- full_join(MetaPhlAn4_tax, tax_df, by = "row.names") %>%
                                    column_to_rownames(var = "Row.names")
        }
      }
    }
  }

  cleanup_tempfiles(c(temp_dir))
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
      original = counts,
      reprocessed = list(
          mOTU3 = mOTU3_counts,
          MetaPhlAn4 = MetaPhlAn4_counts
      )
      ),
      proportions = list(
      original = proportions,
      reprocessed = list(
          mOTU3 = mOTU3_proportions,
          MetaPhlAn4 = MetaPhlAn4_proportions
      )
      ),
      tax = list(
      original = tax,
      reprocessed = list(
          mOTU3 = mOTU3_tax,
          MetaPhlAn4 = MetaPhlAn4_tax
      )
      ),
      scale = scale,
      metadata = metadata
  ))
}