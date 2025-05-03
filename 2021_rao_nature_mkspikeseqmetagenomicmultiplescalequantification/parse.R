parse_2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification <- function(raw = FALSE, originaltax = 'qiime') {
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification")

  # ----- File paths -----
  counts_zip          <- file.path(local, "combined_p1_p2.xlsx.zip")
  supplemental_zip    <- file.path(local, "41586_2021_3241_MOESM4_ESM.xlsx.zip")
  sex_delivery_zip    <- file.path(local, "SI_data3_sex_and_delivery_data.csv.zip")
  diet_data_zip       <- file.path(local, "SI_data2_diet_data.csv.zip")
  meds_data_zip       <- file.path(local, "SI_data1_allMeds_jan2020.xlsx.zip")
  sra_metadata_zip    <- file.path(local, "sra_metadata.csv.zip")
  repro_counts_rds_zip<- file.path(local, "PRJEB36435_dada2_counts.rds.zip")
  repro_tax_zip       <- file.path(local, "PRJEB36435_dada2_taxa.rds.zip")
  
  if (!file.exists(counts_zip)) {
    warning("Counts file not found: ", counts_zip)
  } else {
    tmp_dir <- tempdir()
    unzip(counts_zip, exdir = tmp_dir)
    excel_file <- list.files(tmp_dir, pattern = "combined_p1_p2\\.xlsx$", full.names = TRUE)
    
    sheets <- c("fungi", "bacteria", "archaea", "archaea_side")
    data_list <- lapply(sheets, function(sh) {
      possible_data <- tryCatch(
        read_excel(excel_file, sheet = sh),
        error = function(e) NULL
      )
      return(possible_data)
    })
    names(data_list) <- sheets
    taxonomy_columns <- c(
      "qiime_sklearn", "qiime_confidence",
      "vsearch_usearchglobal", "vsearch_identity",
      "usearch_utax", "usearch_sintax", "usearch_sintax_80%"
    )
    data_fungi    <- data_list$fungi
    data_bacteria <- data_list$bacteria
    if (is.null(data_fungi))    data_fungi <- data.frame()
    if (is.null(data_bacteria)) data_bacteria <- data.frame()
    
    if (nrow(data_bacteria) > 0) {
      # Bacteria
      counts_bacteria <- data_bacteria %>%
        select(OTU_ID, everything()) %>%
        select(-any_of(taxonomy_columns)) 
    }
    
    if (nrow(data_fungi) > 0) {
      # Fungi
      counts_fungi <- data_fungi %>%
        select(OTU_ID, everything()) %>%
        select(-any_of(taxonomy_columns))
    }
}

# ---- Taxa label generator ----
make_taxa_label <- function(df) {
  tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  prefixes  <- c("k", "p", "c", "o", "f", "g")
  if (!all(tax_ranks %in% colnames(df))) {
    stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
  }
  df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
    x[is.na(x) | trimws(x) == ""] <- "unclassified"
    return(x)
  })
  df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
    for (i in length(tax_ranks):1) {
      if (tax_row[i] != "unclassified") {
        return(paste0(prefixes[i], "_", tax_row[i]))
      }
    }
    return("unclassified")  
  })
  return(df)
}

build_taxonomy_table <- function(df, method = c("qiime", "sintax", "sintax_full", "utax", "vsearch")) {
    method <- match.arg(method)

    if (method == "qiime") {
      # pull out the raw qiime_sklearn column
      tax_source <- df %>%
        select(OTU_ID, taxonomy = qiime_sklearn)

      if (any(grepl("^D_\\d+__", tax_source$taxonomy, perl = TRUE))) {
        # 16S style: D_0__..;D_1__..; etc.
        tax_df <- tax_source %>%
          separate(taxonomy,
                  into = c("D0", "D1", "D2", "D3", "D4", "D5", "D6"),
                  sep = ";", fill = "right") %>%
          transmute(
            OTU_ID,
            Kingdom = gsub("^D_0__", "", D0),
            Phylum  = gsub("^D_1__", "", D1),
            Class   = gsub("^D_2__", "", D2),
            Order   = gsub("^D_3__", "", D3),
            Family  = gsub("^D_4__", "", D4),
            Genus   = gsub("^D_5__", "", D5),
            Species = gsub("^D_6__", "", D6)
          )
      } else {
        # ITS style: k__..;p__..; etc.
        tax_df <- tax_source %>%
          separate(taxonomy,
                  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                  sep = ";", fill = "right") %>%
          mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))
      }

      tax_df <- tax_df %>%
        mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), NA, .)))

    } else if (method == "sintax") {
      # usearch_sintax_80% → comma-delimited, no confidences
      tax_df <- df %>%
        select(OTU_ID, taxonomy = `usearch_sintax_80%`) %>%
        separate(taxonomy,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ",", fill = "right") %>%
        mutate(across(Kingdom:Species, ~ gsub("^[a-z]:", "", .)))

    } else if (method == "sintax_full") {
      # usearch_sintax → remove confidence then comma-delimited
      tax_df <- df %>%
        select(OTU_ID, taxonomy = usearch_sintax) %>%
        mutate(taxonomy = gsub("\\([^)]+\\)", "", taxonomy)) %>%
        separate(taxonomy,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ",", fill = "right") %>%
        mutate(across(Kingdom:Species, ~ gsub("^[a-z]:", "", .)))

    } else if (method == "utax") {
      # usearch_utax → semicolon-delimited after last '|refs|'
      tax_df <- df %>%
        select(OTU_ID, taxonomy = usearch_utax) %>%
        mutate(taxonomy = gsub(".*\\|refs\\|", "", taxonomy)) %>%
        separate(taxonomy,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";", fill = "right") %>%
        mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))

    } else if (method == "vsearch") {
      # vsearch_usearchglobal → same as UTAX style
      tax_df <- df %>%
        select(OTU_ID, taxonomy = vsearch_usearchglobal) %>%
        mutate(taxonomy = gsub(".*\\|refs\\|", "", taxonomy)) %>%
        separate(taxonomy,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";", fill = "right") %>%
        mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))
    }

    # final cleanup: blank → NA
    tax_df <- tax_df %>%
      mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), NA, .)))

    return(tax_df)
  }

  taxa_bacteria <- build_taxonomy_table(data_list[["bacteria"]], method = originaltax)
  taxa_fungi    <- build_taxonomy_table(data_list[["fungi"]], method = originaltax)

  taxa_bacteria = make_taxa_label(taxa_bacteria)
  taxa_fungi    = make_taxa_label(taxa_fungi)
  
  tmp_dir <- tempdir()
  unzip(supplemental_zip, exdir = tmp_dir)
  supplemental_excel <- list.files(tmp_dir, pattern = "41586_2021_3241_MOESM4_ESM\\.xlsx$", full.names = TRUE)

  if (length(supplemental_excel) == 0) {
    warning("Supplemental Excel file not found after unzipping.")
  } else {
    file <- supplemental_excel[1]

    read_otu_sheet <- function(sheetname) {
      df <- tryCatch(readxl::read_excel(file, sheet = sheetname), error = function(e) NULL)
      if (is.null(df)) return(NULL)
      if (all(c("qiime2_sklearn_taxonomy", "confidence", "OTU_ID") %in% colnames(df))) {
        tax <- df %>%
          select(OTU_ID, qiime_sklearn = qiime2_sklearn_taxonomy, confidence)
        otu <- df %>%
          select(-qiime2_sklearn_taxonomy, -confidence)
        return(list(otu = otu, tax = build_taxonomy_table(tax, method = "qiime")))
      }
      return(list(otu = df, tax = NULL))
    }

    # Read Scale sheets
    scale_ <- tryCatch(readxl::read_excel(file, sheet = "sTable4"), error = function(e) NULL)
    scale_ITS <- tryCatch(readxl::read_excel(file, sheet = "sTable6"), error = function(e) NULL)
    scale <- bind_rows(Filter(Negate(is.null), list(scale_sheet4, scale_sheet6)))

    # Read OTU sheets
    otu_data <- list(
      originalcounts_bacteria        = read_otu_sheet("sTable8"),
      originalcounts_archaea         = read_otu_sheet("sTable9"),
      originalcounts_fungi           = read_otu_sheet("sTable10"),
      originalcounts_bacteria_phase2 = read_otu_sheet("sTable11"),
      originalcounts_fungi_phase2    = read_otu_sheet("sTable12")
    )

    otu_tables     <- map(otu_data, "otu")
    taxonomy_tables <- map(otu_data, "tax")

    # ---- Load Sheet3 and extract mock_metadata ----
    sheet3 <- read_excel(file, sheet = "sTable3", col_names = FALSE)

    mock_metadata <- sheet3 %>%
      slice(3:which(.data[[1]] == "3b. OTU table of bacterial mock communities") - 1) %>%
      {
        colnames(.) <- as.character(sheet3[2, 1:ncol(.)])
        .
      } %>%
      rename_with(~ make.names(., unique = TRUE)) %>%
      mutate(across(everything(), readr::parse_guess))

    # ---- Function to read each OTU table from its labeled section ----
    read_mock_table <- function(label, n_max) {
      start <- which(sheet3[[1]] == label) + 2
      tbl <- read_excel(file, sheet = "sTable3", skip = start, n_max = n_max)
      tbl <- tbl[1:(n_max - 1), ]
      return(tbl)
    }

    # ---- Load the bacterial, fungal, and fecal test OTU tables ----
    mock_otu_tables <- list(
      originalcounts_mockbacteria = read_mock_table("3b. OTU table of bacterial mock communities", 13),
      originalcounts_mockfungi    = read_mock_table("3c. OTU table of fungal mock communities", 12),
      originalcounts_mockfecal  = read_mock_table("3d. OTU table of test fecal samples", 43)
    )

    testfecal_labeled <- make_taxa_label(mock_otu_tables$originalcounts_mockfecal)
    originaltax_mockfecal <- testfecal_labeled %>%
      select(Kingdom, Phylum, Class, Order, Family, Genus, Taxa) %>%
      as.data.frame()

    mock_otu_tables$originalcounts_mockfecal <- testfecal_labeled %>%
      column_to_rownames("Taxa") %>%
      select(-c(Kingdom, Phylum, Class, Order, Family, Genus)) %>%
      t() %>%
      as.data.frame()

    library(tidyverse)

    # ---- Taxonomy table constructor ----
    make_mock_tax_table <- function(df) {
      df <- df %>%
        filter(!taxonomy %in% c("Observed_total", "Expected_total")) %>%
        mutate(
          Genus   = sub("\\..*", "", taxonomy),
          Species = sub(".*?\\.\\s*", "", taxonomy),
          Kingdom = NA, Phylum = NA, Class = NA, Order = NA, Family = NA
        ) %>%
        mutate(across(everything(), ~ ifelse(. == "", NA, .))) %>%
        mutate(Taxa = paste0("g_", Genus, "_", Species)) %>%
        select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Taxa)
      
      return(df)
    }

    # ---- Process one mock OTU table (taxonomy + scale extraction + OTU formatting) ----
    process_mock_table <- function(df) {
      # Separate scale rows
      scale_df <- df %>%
        filter(taxonomy %in% c("Observed_total", "Expected_total")) %>%
        column_to_rownames("taxonomy") %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("Sample") %>%
        rename(Observed_total = Observed_total, Expected_total = Expected_total)

      # Process taxonomy table
      tax_df <- make_mock_tax_table(df)

      # Subset OTU table (excluding scale rows)
      otu_df <- df %>%
        filter(!taxonomy %in% c("Observed_total", "Expected_total")) %>%
        mutate(Taxa = paste0("g_", sub("\\..*", "", taxonomy), "_", sub(".*?\\.\\s*", "", taxonomy))) %>%
        select(-taxonomy) %>%
        column_to_rownames("Taxa")

      return(list(otu = otu_df, tax = tax_df, scale = scale_df))
    }

    # ---- Apply to both mockbacteria and mockfungi ----
    mock_bacteria_result <- process_mock_table(mock_otu_tables$originalcounts_mockbacteria)
    mock_fungi_result    <- process_mock_table(mock_otu_tables$originalcounts_mockfungi)
    mock_otu_tables$originalcounts_mockbacteria <- as.data.frame(t(mock_bacteria_result$otu))
    mock_otu_tables$originalcounts_mockfungi    <- as.data.frame(t(mock_fungi_result$otu))
    originaltax_mockbacteria <- mock_bacteria_result$tax
    originaltax_mockfungi    <- mock_fungi_result$tax
    scale_mockbacteria <- mock_bacteria_result$scale
    scale_mockfungi    <- mock_fungi_result$scale

  }

  make_proportions <- function(count_df) {
    if (is.null(count_df) || nrow(count_df) == 0) return(NULL)

    count_df_num <- count_df %>%
      mutate(across(everything(), as.numeric))

    row_sums <- rowSums(count_df_num, na.rm = TRUE)
    prop_df <- sweep(count_df_num, 1, row_sums, "/")
    prop_df[is.na(prop_df)] <- 0
    return(prop_df)
  }

  proportions_bacteria                = make_proportions(counts_bacteria)
  proportions_fungi                   = make_proportions(counts_fungi)
  originalproportions_bacteria        = make_proportions(originalcounts_bacteria)
  originalproportions_archaea         = make_proportions(originalcounts_archaea)
  originalproportions_fungi           = make_proportions(originalcounts_fungi)
  originalproportions_bacteria_phase2 = make_proportions(originalcounts_bacteria_phase2)
  originalproportions_fungi_phase2    = make_proportions(originalcounts_fungi_phase2)
  originalproportions_mockbacteria    = make_proportions(mock_otu_tables$originalcounts_mockbacteria)
  originalproportions_mockfungi       = make_proportions(mock_otu_tables$originalcounts_mockfungi)
  originalproportions_mockfecal       = make_proportions(mock_otu_tables$originalcounts_mockfecal)

  metadata <- NULL
  
  safe_read_zip <- function(zip_path, is_xlsx = FALSE, sheet = 1) {
    if (!file.exists(zip_path)) return(NULL)
    zfiles <- unzip(zip_path, list = TRUE)$Name
    if (length(zfiles) < 1) return(NULL)
    file_to_read <- zfiles[1]
    td <- tempdir()
    unzip(zip_path, files = file_to_read, exdir = td, overwrite = TRUE)
    fpath <- file.path(td, file_to_read)
    if (!file.exists(fpath)) return(NULL)
    df <- tryCatch(
      {
        if (is_xlsx) {
          read_excel(fpath, sheet = sheet) %>%
            as.data.frame()
        } else {
          read.csv(fpath, stringsAsFactors = FALSE) %>%
            as.data.frame()
        }
      },
      error = function(e) NULL
    )
    return(df)
  }
  
  # Read existing metadata files
  sex_delivery_data <- safe_read_zip(sex_delivery_zip, is_xlsx = FALSE)
  diet_data         <- safe_read_zip(diet_data_zip, is_xlsx = FALSE)
  meds_data         <- safe_read_zip(meds_data_zip, is_xlsx = TRUE, sheet = 1)
  sra_metadata      <- safe_read_zip(sra_metadata_zip, is_xlsx = FALSE)
  
  if (!is.null(sex_delivery_data)) {
    if ("baby_id" %in% colnames(sex_delivery_data)) {
      sex_delivery_data <- sex_delivery_data %>%
        rename(id = baby_id)
    }
    metadata <- sex_delivery_data
    if (!is.null(diet_data)) {
      metadata <- metadata %>%
        full_join(diet_data, by = "id")
    }
    if (!is.null(meds_data)) {
      metadata <- metadata %>%
        full_join(meds_data, by = "id")
    }
  }
  
  # Note: The sample_name column in sra_metadata (e.g., "Ca-10-1-R1-low_S70_L001_R1") may not directly match the 'id' column in the existing metadata.
  # may need to clean and standardize the sample_name and id columns to ensure proper merging (e.g., by extracting a common identifier or removing suffixes).
  
  # Merge SRA metadata using sample_name
  if (!is.null(sra_metadata) && !is.null(metadata)) {
    # Attempt to merge using sample_name and id
    metadata <- metadata %>%
      full_join(sra_metadata, by = c("id" = "sample_name"))
    metadata <- as.data.frame(metadata)
  }


  # NEED TO SEPARATE THE DIFFERENT SAMPLES INTO SPECIFIC DATASETS

    # ----- Reprocessed counts from RDS ZIP -----
  temp_rds <- tempfile(fileext = ".rds")
  unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

  rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
  counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

  # ----- Taxonomy reprocessed -----
  temp_tax <- tempfile(fileext = ".rds")
  unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

  tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
  if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
  tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

  
  # ----- Convert sequences to lowest rank taxonomy found and update key -----
  make_taxa_label <- function(df) {
      tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
      prefixes  <- c("k", "p", "c", "o", "f", "g")
      if (!all(tax_ranks %in% colnames(df))) {
      stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
      }
      df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
      x[is.na(x) | trimws(x) == ""] <- "unclassified"
      return(x)
      })
      df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
      for (i in length(tax_ranks):1) {
          if (tax_row[i] != "unclassified") {
          return(paste0(prefixes[i], "_", tax_row[i]))
          }
      }
      return("unclassified")  
      })
      return(df)
  }
  tax_reprocessed = make_taxa_label(tax_reprocessed)

  # ----- Convert accessions to sample IDs / Sequences to Taxa -----
  # accessions to sampleIDs is study specific: IF NEED BE

  # taxa
  if (!raw) {
      matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
      colnames(counts_reprocessed) <- matched_taxa
      counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
  }

  # proportions reprocessed
  proportions_reprocessed = counts_reprocessed
  proportions_reprocessed[-1] <- lapply(
      counts_reprocessed[-1],
      function(col) col / sum(col)
  )
  
  # -------------------------------------------------------------
  return(list(
    counts = list(
      original = list(
        combinedphases = list(
          bacteria = counts_bacteria,
          archaea = NA,
          fungi = counts_fungi
        ),
        phase1 = list(
          bacteria = originalcounts_bacteria,
          archaea = originalcounts_archaea,
          fungi = originalcounts_fungi
        ),
        phase2 = list(
          bacteria = originalcounts_bacteria_phase2,
          archaea = NA,
          fungi = originalcounts_fungi_phase2
        ),
        mock = list(
          bacteria = mock_otu_tables$originalcounts_mockbacteria,
          archaea = NA,
          fungi = mock_otu_tables$originalcounts_mockfungi,
          fecal = mock_otu_tables$originalcounts_mockfecal
        )
      ),
      reprocessed = list(

      )
    ),
    proportions = list(
      original = list(
        combinedphases = list(
          bacteria = proportions_bacteria,
          archaea = NA,
          fungi = proportions_fungi
        ),
        phase1 = list(
          bacteria = originalproportions_bacteria,
          archaea = originalproportions_archaea,
          fungi = originalproportions_fungi
        ),
        phase2 = list(
          bacteria = originalproportions_bacteria_phase2,
          archaea = NA,
          fungi = originalproportions_fungi_phase2
        ),
        mock = list(
          bacteria = mock_otu_tables$originalproportions_mockbacteria,
          archaea = NA,
          fungi = mock_otu_tables$originalproportions_mockfungi,
          fecal = mock_otu_tables$originalproportions_mockfecal
        )
      ),
      reprocessed = list(

      )
    ),
    tax = list(
      original = list(
        combinedphases = list(
          bacteria = proportions_bacteria,
          archaea = NA,
          fungi = proportions_fungi
        ),
        phase1 = list(
          bacteria = originalproportions_bacteria,
          archaea = originalproportions_archaea,
          fungi = originalproportions_fungi
        ),
        phase2 = list(
          bacteria = originalproportions_bacteria_phase2,
          archaea = NA,
          fungi = originalproportions_fungi_phase2
        ),
        mock = list(
          bacteria = mock_otu_tables$originalproportions_mockbacteria,
          archaea = NA,
          fungi = mock_otu_tables$originalproportions_mockfungi,
          fecal = mock_otu_tables$originalproportions_mockfecal
        )
      ),
      reprocessed = list(

      )
    ),
    scale = scale,
    metadata = metadata
  ))
}