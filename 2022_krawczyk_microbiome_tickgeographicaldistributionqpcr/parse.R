parse_2022_krawczyk_microbiome_tickgeographicaldistributionqpcr <- function(raw=FALSE) {
  required_pkgs <- c("tibble", "tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)
  library(readxl)

  # ----- Local base directory -----
  local <- file.path("2022_krawczyk_microbiome_tickgeographicaldistributionqpcr")

  # ----- File paths -----
  metadata_two_zip     <- file.path(local, "40168_2022_1276_MOESM1_ESM.zip")
  metadata_zip         <- file.path(local, "SraRunTable (30).csv.zip")
  counts_zip           <- file.path(local, "40168_2022_1276_MOESM3_ESM.zip")
  scale_zip            <- file.path(local, "scale_qpcr.zip")
  repro_counts_rds_zip <- file.path(local, "PRJNA813158_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA813158_dada2_taxa.rds.zip")

  # ----- Scale -----
  scale_file <- unzip(scale_zip, list = TRUE)$Name[1]
  scale_con  <- unz(scale_zip, scale_file)
  scaledata  <- read.csv(scale_con) %>% as.data.frame() 
  cleaned <- as.numeric(scaledata$`X16S.rRNA.content.in.ng.µL`)  # Will coerce character "no DNA left" to NA
  scale = data.frame(Sample_name = scaledata$`Sample.ID`, `16S rRNA content in ng/µL` = cleaned)
  
  # ----- Metadata -----
  meta_file <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con  <- unz(metadata_zip, meta_file)
  metadata  = read.csv(meta_con) %>% as.data.frame() %>% rename(Sample_name = `Sample.ID`)
  meta_two_file <- unzip(metadata_two_zip, list = TRUE)$Name[1]
  extracted_xlsx <- unzip(metadata_two_zip, files = meta_two_file, exdir = tempdir(), overwrite = TRUE)[1]
  metadata_two <- read_excel(extracted_xlsx, sheet = 2) %>% rename(Sample_name = `Sample ID`)
  metadata <- left_join(metadata, metadata_two, by = "Sample_name") %>% 
              rename(Accession = Run)

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
  make_taxa_label <- function(df) {
    tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    prefixes  <- c("k", "p", "c", "o", "f", "g")
    if (!all(tax_ranks %in% colnames(df))) {
        stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
    }
    df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
        x[is.na(x) | trimws(x) == ""] <- "unclassified"
        x
    })
    df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
        if (tax_row["Genus"] != "unclassified") {
        return(paste0("g_", tax_row["Genus"]))
        }
        for (i in (length(tax_ranks)-1):1) {  # skip Genus
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

  # -------- Counts --------
  if (file.exists(counts_zip)) {
    counts_file <- unzip(counts_zip, list = TRUE)$Name[1]
    extracted_xlsx <- unzip(counts_zip, files = counts_file, exdir = tempdir(), overwrite = TRUE)[1]
    counts_original = read_excel(extracted_xlsx, sheet = 1)
    columns_to_drop <- c("Name","Taxonomy", "Combined Abundance",	"Min",	"Max",	"Mean",	"Median",	"Std")

    # Taxa
    taxonomy_cols <- c("Sequence","Taxonomy")
    tax <- counts_original[, names(counts_original) %in% taxonomy_cols]
    tax <- tax %>%
      separate(
        col = Taxonomy,
        into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
        sep = ",\\s*",
        remove = FALSE
      ) %>%
      mutate(across(Kingdom:Genus, ~ sub("D_\\d+__", "", .)))
    tax = make_taxa_label(tax)

    # Counts
    rownames(tax) <- tax$Sequence
    counts_original <- counts_original %>% column_to_rownames("Sequence")
    counts_original = counts_original[, !(names(counts_original) %in% columns_to_drop)]
    counts_original = as.data.frame(t(counts_original))

    if (!raw) {
      matched_taxa <- tax$Taxa[match(colnames(counts_original), rownames(tax))]
      colnames(counts_original) <- matched_taxa
      counts_original <- as.data.frame(t(rowsum(t(counts_original), group = colnames(counts_original))))
    }
  } else {
    counts_original = NA
  }

  proportions_original = counts_original
  proportions_original[-1] <- lapply(
    counts_original[-1],
    function(col) col / sum(col)
  )

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

  # DELETE LATER #####################################
  maxwellreprocessedpreviously = file.path(local, "Krawczyk_2022_16S.csv.zip")
  counts_original = read_zipped_table(maxwellreprocessedpreviously, row.names = NULL) %>% as.data.frame()
  proportions_original <- sweep(counts_original, MARGIN = 1,STATS  = rowSums(counts_original), FUN = "/")
  ####################################################
  if (!raw) {
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

  # ----- Return all -----
  return(list(
    counts      = list(
      original = counts_original,
      reprocessed = counts_reprocessed
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = proportions_reprocessed
    ),
    tax         = list(
      original = tax,
      reprocessed = tax_reprocessed
    ),
    scale       = scale,
    metadata    = metadata
  ))
}
