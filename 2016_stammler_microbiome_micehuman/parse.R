parse_2016_stammler_microbiome_micehuman <- function(raw = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2016_stammler_microbiome_micehuman")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJEB11953_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJEB11953_dada2_taxa.rds.zip")
    scale_16s_zip        <- file.path(local, "Stammler2016_scale.csv.zip")
    counts_16s_zip       <- file.path(local, "Stammler_2016_16S.csv.zip")
    metadata_16s_zip     <- file.path(local, "Stammler_2016_metadata.csv.zip")
    metadata_SRA_zip     <- file.path(local, "SraRunTable (38).csv.zip")

    read_zipped_csv <- function(zip_path) {
      if (file.exists(zip_path)) {
          csv_file <- unzip(zip_path, list = TRUE)$Name[1]
          read.csv(unz(zip_path, csv_file), row.names = 1, check.names = FALSE)
      } else {
          warning(paste("File not found:", zip_path))
          return(NA)
      }
    }

    # ----- Initialize everything as NA -----
    counts_original_mice <- NA
    proportions_original_mice <- NA
    tax_original_mice <- NA
    counts_original_human <- NA
    proportions_original_human <- NA
    tax_original_human <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA

    # ------ original counts ------
    counts_original_mice <- read_zipped_csv(counts_16s_zip)

    if (!is.na(counts_original_mice)[1]) {
      original_taxa <- colnames(counts_original_mice)

      # Create taxa mapping data frame
      tax_original_mice <- data.frame(
        Taxa = original_taxa,
        stringsAsFactors = FALSE
      )

      # ------ proportions from counts ------
      proportions_original_mice <- sweep(
        counts_original_mice,
        MARGIN = 1,                          
        STATS = rowSums(counts_original_mice, na.rm = TRUE),
        FUN = "/"
      )
    } else {
      proportions_original_mice <- NA
      tax_original_mice <- NA
    }

    # ---- scale and metadata -----
    scale     <- read_zipped_csv(scale_16s_zip)
    metadata  <- read_zipped_csv(metadata_16s_zip)
    sra       <- read_zipped_csv(metadata_SRA_zip)

    sra <- sra %>%
      mutate(Sample_name = sub(".*(MID.*)", "\\1", Sample_name))

    metadata  =  metadata %>%
      full_join(sra, by = c("SampleID" = "Sample_name")) %>%
      rename(Accession = Run)
    scale     = scale %>%
      left_join(metadata %>% select(Sample_name, Accession), by = c("SampleID" = "Sample_name"))

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

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original_mice,
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = proportions_original_mice,
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = tax_original_mice,
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))
}



