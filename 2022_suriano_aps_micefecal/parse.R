parse_2022_suriano_aps_micefecal <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "readr")
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
    library(readr)

    # ----- local path ------------------------
    local                   <- file.path("2022_suriano_aps_micefecal")

    # ----- file paths -------------------------
    metadata_zip            <- file.path(local,"SraRunTable.csv.zip")
    scale_zip               <- file.path(local,"Supplementary Table 6.xlsx.zip")
    repro_counts_rds_zip    <- file.path(local,"PRJEB53668_dada2_counts.rds.zip")
    repro_tax_zip           <- file.path(local,"PRJEB53668_dada2_taxa.rds.zip")

    # ----- metadata ---------------------------
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]  # list file inside zip
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
    metadata$Sample_name <- gsub("^(.*?_D)(\\d{2})\\.(\\d+)$", "D\\2-\\3", metadata$Sample_name)
    metadata = metadata %>% rename(Accession = Run)
    
    # ----- scale ------------------------------
    scale_xlsx <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_xlsx, exdir = tempdir(), overwrite = TRUE)
    
    scale <- read_excel(scale_path, sheet = "Microbial loads")
    scale <- scale %>%
        select(Sample, `Cells.g.of.fecal.sample`) %>%
        rename(Sample_name = Sample) %>%
        mutate(log2_FC_cells_per_g = ifelse(`Cells.g.of.fecal.sample` > 0, log2(`Cells.g.of.fecal.sample`), NA)) %>%
        mutate(log10_FC_cells_per_g = ifelse(`Cells.g.of.fecal.sample` > 0, log10(`Cells.g.of.fecal.sample`), NA))
    
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
    if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
    }

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
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }

    return(list(
        counts = list(
            original = NA,
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = NA,
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = NA,
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))
}