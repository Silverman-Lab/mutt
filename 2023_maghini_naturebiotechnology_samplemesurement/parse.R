parse_2023_maghini_naturebio_metagenomic <- function(raw = FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "readr", "taxizedb")
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
  library(taxizedb)

  # ----- Local base directory -----
  local <- file.path("2023_maghini_naturebiotechnology_samplemesurement")

  # ----- File paths -----
  motus_zip      <- file.path(local, "PRJNA940499_motus_merged.tsv.zip")
  metaphlan4_zip <- file.path(local, "PRJNA940499_MetaPhlAn_merged.tsv.zip")
  scale_zip      <- file.path(local, "Maghini2023_scale.csv.zip")
  metadata_zip   <- file.path(local, "Maghini_2023_metadata.csv.zip")
  original_zip   <- file.path(local, "Maghini_2023_shotgunmetagenomics.csv.zip")

  read_zipped_table <- function(zip_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE) {
        if (file.exists(zip_path)) {
        inner_file <- unzip(zip_path, list = TRUE)$Name[1]
        con <- unz(zip_path, inner_file)
        read.table(con, sep = sep, header = header, row.names = row.names, check.names = check.names, stringsAsFactors = FALSE)
        } else {
        warning(paste("File not found:", zip_path))
        return(NA)
        }
    }

  # ----- Convert sequences to lowest rank taxonomy found and update key -----
  make_taxa_label <- function(df) {
          tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
          prefixes  <- c("k", "p", "c", "o", "f", "g", "s")
          if (!all(tax_ranks %in% colnames(df))) {
              stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
          }
          df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
              x[is.na(x) | trimws(x) == ""] <- "unclassified"
              x
          })
          df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
              if (tax_row["Species"] != "unclassified") {
              return(paste0("s_", tax_row["Species"]))
              }
              for (i in (length(tax_ranks)-1):1) {  
              if (tax_row[i] != "unclassified") {
                  return(paste0("uc_", prefixes[i], "_", tax_row[i]))
              }
              }
              return("unclassified")
          })
          return(df)
      }

  fill_na_zero_numeric <- function(x) {
      if (is.data.frame(x)) {
          x[] <- lapply(x, function(y) if (is.numeric(y)) replace(y, is.na(y), 0) else y)
      } else if (is.matrix(x) && is.numeric(x)) {
          x[is.na(x)] <- 0
      } else if (is.list(x)) {
          x <- lapply(x, fill_na_zero_numeric)
      }
      x
  }

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

  # ----- Read taxonomic counts -----

  ## UNCOMMENT WHEN FIXED 

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
  #     matched_taxa <- tax$Taxa[match(colnames(counts_original), rownames(tax))]
  #     colnames(counts_original) <- matched_taxa
  #     counts_original <- as.data.frame(t(rowsum(t(counts_original), group = colnames(counts_original))))
  # }

  # # ------- Proportions ----------
  # proportions_original <- map(counts_original, function(df) {
  #   rsums <- rowSums(df)
  #   prop  <- sweep(df, 1, rsums, FUN = "/")
  #   prop[is.nan(prop)] <- 0
  #   return(as_tibble(prop, rownames = "Taxon"))
  # })

  # ----- Read and filter qPCR data -----
  plate_list <- list.files(path = base_path, pattern = "qPCR_plate", full.names = TRUE)
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
    )

  # ----- Extract and deduplicate metadata -----
  metadata <- qPCR %>%
    dplyr::select(ID, Donor, Condition, Replicate) %>%
    distinct()

  # ----- SRA metadata and merge -----
  metadatafromsra_path <- paste0(base_path, "/SraRunTable (3).csv")
  metadatafromsra <- read.csv(metadatafromsra_path, check.names = FALSE) %>%
    mutate(Library.Name = gsub("(_DNA)", "", Library.Name)) %>% rename(Accession = Run)

  scale$ID <- gsub("-", "_", scale$ID)
  metadata$ID <- gsub("-", "_", metadata$ID)

  metadata <- merge(
    metadatafromsra, metadata,
    by.x = "Library.Name",
    by.y = "ID",
    all = TRUE
  )

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
        proportions <- apply(df, 2, function(col) col / sum(col))
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
        proportions <- apply(df, 2, function(col) col / sum(col))
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

  # Can delete later:
  counts_original = read_zipped_table(original_zip)
  proportions_original <- map(counts_original, function(df) {
    rsums <- rowSums(df)
    prop  <- sweep(df, 1, rsums, FUN = "/")
    prop[is.nan(prop)] <- 0
    return(as_tibble(prop, rownames = "Taxon"))
  })
  
  
  # Dont delete:
  if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      #mOTU3_counts = fill_na_zero_numeric(mOTU3_counts)
      proportions_original = fill_na_zero_numeric(proportions_original)
      #MetaPhlAn4_counts = fill_na_zero_numeric(MetaPhlAn4_counts)
      #mOTU3_proportions = fill_na_zero_numeric(mOTU3_proportions)
      #MetaPhlAn4_proportions = fill_na_zero_numeric(MetaPhlAn4_proportions)
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


