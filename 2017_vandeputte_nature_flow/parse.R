parse_2017_vandeputte_nature_flow <- function(raw = FALSE) {
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
    local <- file.path("2017_vandeputte_nature_flow")

    # ----- File paths -----
    metadata_zip         <- file.path(local, "Vandeputte_2017_metadata.csv.zip")
    metadata_two_zip     <- file.path(local, "cellcountstotal.csv.zip")
    orig_counts_zip      <- file.path(local, "OTU_nochim.zip")
    orig_tax_rdp_zip     <- file.path(local, "otu_taxonomy_rdp.csv.zip")
    orig_tax_silva_zip   <- file.path(local, "otu_taxonomy_silva.csv.zip")
    orig_prop_zip        <- NA
    repro_counts_rds_zip <- file.path(local, "PRJEB21504_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJEB21504_dada2_taxa.rds.zip")


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
      if (missing(x)) return(NULL)
      if (is.data.frame(x)) {
          x[] <- lapply(x, function(y) if (is.numeric(y)) replace(y, is.na(y), 0) else y)
      } else if (is.matrix(x) && is.numeric(x)) {
          x[is.na(x)] <- 0
      } else if (is.list(x)) {
          x <- lapply(x, fill_na_zero_numeric)
      }
      x
    }
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

    # ----- Metadata and Scale -----
    # Read scale
    count_csv    <- unzip(metadata_two_zip, list = TRUE)$Name[1]
    count_con    <- unz(metadata_two_zip, count_csv)
    metadata_two <- read.csv(count_con) %>% as.data.frame()

    # Read metadata
    meta_csv     <- unzip(metadata_zip, list = TRUE)$Name[1]
    meta_con     <- unz(metadata_zip, meta_csv)
    metadata_df  <- read.csv(meta_con) %>% as.data.frame()

    # Join and select
    df <- full_join(
    metadata_df,
    metadata_two,
    by = c("Sample" = "Individual")
    )
    metadata <- df %>%
    select(
      Accession = sampleID,
      Sample,
      Cohort = Cohort.x,
      Day,
      `Health status` = Health.status,
      Enterotype = Enterotype.x
    )
    scale <- df %>%
    select(
      Sample,
      Day,
      Accession = sampleID,
      `Average cell count (per gram of fresh feces)` = Average.cell.count..per.gram.of.fresh.feces.,
      `STDEV cell count (per gram of fresh feces)` = STDEV.cell.count..per.gram.of.fresh.feces.,
      `Average cell count (per gram of frozen feces)` = Average.cell.count..per.gram.of.frozen.feces..y,
      `STDEV cell count (per gram of frozen feces)` = STDEV.cell.count..per.gram.of.frozen.feces..y
    )

    # ----- Original counts from CSV.zip -----
    if (file.exists(orig_counts_zip)) {
    orig_csv <- unzip(orig_counts_zip, list = TRUE)$Name[1]
    orig_con <- unz(orig_counts_zip, orig_csv)
    orig_mat <- read.csv(orig_con, row.names = 1, check.names = FALSE)
    counts_original <- as.data.frame(orig_mat)
    counts_original$Sequence <- rownames(counts_original)
    counts_original <- counts_original[, c("Sequence", setdiff(names(counts_original), "Sequence"))]
    rownames(counts_original) <- paste0("Taxon_", seq_len(nrow(counts_original)))
    } else {
    counts_original <- NA
    }

    if (file.exists(orig_counts_zip)) {
    proportions_original = counts_original
    proportions_original[-1] <- lapply(
      counts_original[-1],
      function(col) col / sum(col)
    )
    } else if (file.exists(orig_prop_zip)) {
    prop_csv <- unzip(orig_prop_zip, list = TRUE)$Name[1]
    prop_con <- unz(orig_prop_zip, prop_csv)
    proportions_original = read.csv(prop_con, row.names = 1, check.names = FALSE) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sequence") %>%
      dplyr::select(Sequence, everything())
    } else {
    proportions_original <- NA
    }

    # --- Original taxonomies ---
    read_taxonomy_zip <- function(zip_path) {
    if (!file.exists(zip_path)) return(NA)
    zip_contents <- unzip(zip_path, list = TRUE)
    csv_name <- zip_contents$Name[1]
    con <- unz(zip_path, csv_name)
    tax_df <- read.csv(con, row.names = 1, check.names = FALSE)
    tax_df$Sequence <- rownames(tax_df)
    tax_df <- tax_df[, c("Sequence", setdiff(names(tax_df), "Sequence"))]
    }

    tax_original_rdp   <- read_taxonomy_zip(orig_tax_rdp_zip)
    tax_original_silva <- read_taxonomy_zip(orig_tax_silva_zip)

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

    if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_original = fill_na_zero_numeric(proportions_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }

    # ----- Return all -----
    return(list(
      counts      = list(
        original    = counts_original,
        reprocessed = counts_reprocessed
      ),
      proportions = list(
        original    = proportions_original,
        reprocessed = proportions_reprocessed
      ),
      tax         = list(
        original_rdp    = tax_original_rdp,
        original_silva  = tax_original_silva,
        reprocessed = tax_reprocessed
      ),
      scale       = scale,
      metadata    = metadata
    ))
}
