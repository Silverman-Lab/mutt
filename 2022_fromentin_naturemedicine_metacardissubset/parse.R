parse_2022_fromentin_naturemedicine_metacardissubset <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("tidyverse", "readxl", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
      stop(
          "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
          ". Please install them before running this function."
      )
  }
  # Load needed libraries
  library(tidyverse)
  library(readxl)
  library(readr)

  # -------- local path -----------
  local                       <- file.path("2021_forslund_nature_metacardis")

  # -------- file paths -----------
  metamatmetformin_zip        <- file.path(local, "metaMatMetformin (1).RData.zip")
  metadata_zip                <- file.path(local, "41591_2022_1688_MOESM3_ESM.xlsx.zip")
  sra_zips                    <- c(
    file.path(local, "SraRunTable (36).csv.zip")
  )
  repro_motus_zips            <- c(
    file.path(local, "PRJEB46098_motus_merged.tsv.zip")  
  )

  repro_metaphlan_zips        <- c(
    file.path(local, "PRJEB46098_MetaPhlAn_merged.tsv.zip")
  )

  if (file.exists(metadata_zip)) {
    temp_dir <- tempfile("metadata_unzip_")
    dir.create(temp_dir)
    unzip(metadata_zip, exdir = temp_dir)

    excel_files <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)
    
    if (length(excel_files) > 0) {
      metadata_file <- excel_files[1]
      metadata_sheets <- c("ST9","ST10", "ST11", "ST12", "ST13", "ST14")
      
      metadata_list <- list()
      scale <- NULL
      st10_selected <- NULL
      st14_df <- NULL

      for (sheet in metadata_sheets) {
        df_raw <- read_excel(metadata_file, sheet = sheet, col_names = FALSE)
        if (nrow(df_raw) < 3) next
        df_name <- as.character(df_raw[[1, 1]])
        colnames(df_raw) <- as.character(unlist(df_raw[2, ]))
        df <- df_raw[-c(1, 2), , drop = FALSE]
        df <- as.data.frame(df)
        if ("ID" %in% colnames(df)) {
          df <- df %>%
            rename(Sample = ID) %>%
            mutate(
              Sample = if_else(
                grepl("^x", Sample, ignore.case = FALSE),
                paste0("M0", Sample),
                Sample
              )
            )
        }
        if (sheet == "ST10") {
          if (all(c("Sample", "Microbial load") %in% colnames(df))) {
            scale <- df %>%
              select(Sample, MicrobialLoad = `Microbial load`) %>%
              as.data.frame(stringsAsFactors = FALSE)
          }
          if (all(c("Status", "MGS count", "Gene count", "Microbial load") %in% colnames(df))) {
            st10_selected <- df %>%
              select(Sample, Status, `MGS count`, `Gene count`, `Microbial load`)
          }
        } else if (sheet == "ST14") {
          st14_df <- df
        } else if (sheet == "ST9") {
          st9_df <- df
        } else {
          df <- df %>% select(-any_of(c("Status", "MGS count", "Gene count", "Microbial load")))
          metadata_list[[df_name]] <- df
        }
      }

      if (!is.null(st10_selected) && !is.null(st14_df) && !is.null(st9_df)) {
        metadata_df <- list(st14_df, st10_selected, st9_df) %>%
                        reduce(full_join, by = "Sample")
      } else {
        metadata_df <- NULL
      }
    } else {
      warning("No .xlsx file found in ", metadata_zip)
    }

    unlink(temp_dir, recursive = TRUE)
  } else {
    warning("Metadata zip file not found: ", metadata_zip)
  }

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

  metadata <- full_join(
    combined_sra_df,
    metadata_df,
    by = c("Submitter_Id" = "Sample") 
  )

  # --- original counts, proportions, tax ----

  # NEED TO FINALIZE THIS

  # ---- initialize dataframes ----
  counts = NA
  proportions = NA
  tax = NA

  # reduced_feature_zip <- paste0(local, "reduced_feature (1).RData.zip")
  # if (file.exists(reduced_feature_zip)) {
    
  #   temp_dir <- tempdir()
  #   unzip(reduced_feature_zip, exdir = temp_dir)
  #   rdata_files <- list.files(temp_dir, pattern = "\\.RData$", full.names = TRUE)
    
  #   for (rdata_file in rdata_files) {
  #     env <- new.env()
  #     load(rdata_file, envir = env)

  #     if (exists("reduced_feature", envir = env)) {
  #       out$counts <- env$reduced_feature
  #     } else {
  #       warning("reduced_feature not found in ", rdata_file)
  #     }
  #   }
  # } else {
  #   warning("Reduced_feature zip file not found: ", reduced_feature_zip)
  # }

  # metamatmetformin_zip <- paste0(local, "metaMatMetformin (1).RData.zip")
  # if (file.exists(metamatmetformin_zip)) {
  #   message("Extracting and loading metamatmetformin data from ", metamatmetformin_zip)
    
  #   temp_dir <- tempdir()
  #   unzip(metamatmetformin_zip, exdir = temp_dir)
  #   rdata_files <- list.files(temp_dir, pattern = "\\.RData$", full.names = TRUE)
    
  #   for (rdata_file in rdata_files) {
  #     env <- new.env()
  #     load(rdata_file, envir = env)

  #     if (exists("metamatmetformin", envir = env)) {
  #       out$metadata_small <- env$metamatmetformin
  #     } else {
  #       warning("metamatmetformin not found in ", rdata_file)
  #     }
  #   }
  # } else {
  #   warning("Metamatmetformin zip file not found: ", metamatmetformin_zip)
  # 

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
          aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
          df = aligned$reprocessed
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
          aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
          df = aligned$reprocessed
        }
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



  # ----- Return -----
  return(list(
      counts = list(
        original = list(
            microbiome = if (!raw) {fill_na_zero_numeric(metadata_list$microbiome)} else {metadata_list$microbiome},
            GMMandKEGGabundances = metadata_list$GMMandKEGGabundances,
            serummetabolites = NA,
            urinemetabolites = metadata_list$urinemetabolitesdata
        ),
        reprocessed = list(
            mOTU3 = if (!raw) {fill_na_zero_numeric(mOTU3_counts)} else {mOTU3_counts},
            MetaPhlAn4 = if (!raw) {fill_na_zero_numeric(MetaPhlAn4_counts)} else {MetaPhlAn4_counts}
        )
      ),
      proportions = list(
        original = list(
            microbiome = if (!raw) {fill_na_zero_numeric(proportions_list$microbiome)} else {proportions_list$microbiome},
            GMMandKEGGabundances = proportions_list$GMMandKEGGabundances,
            serummetabolites = metadata_list$logtransformedserummetabolites,
            urinemetabolites = proportions_list$urinemetabolites
        ),
        reprocessed = list(
            mOTU3 = if (!raw) {fill_na_zero_numeric(mOTU3_proportions)} else {mOTU3_proportions},
            MetaPhlAn4 = if (!raw) {fill_na_zero_numeric(MetaPhlAn4_proportions)} else {MetaPhlAn4_proportions}
        )
      ),
      tax = list(
        original = list(
            microbiome = tax_list$microbiome,
            GMMandKEGGabundances = tax_list$GMMandKEGGabundances,
            serummetabolites = tax_list$logtransformedserummetabolites,
            urinemetabolites = tax_list$urinemetabolites
        ),
        reprocessed = list(
            mOTU3 = mOTU3_tax,
            MetaPhlAn4 = MetaPhlAn4_tax
        )
      ),
      scale = scale,
      metadata = metadata
  ))
}