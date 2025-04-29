parse_2020_vierasilva_nature_BMIS <- function() {
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
    local                       <- file.path("2020_vierasilva_nature_BMIS")

    # -------- file paths -----------
    original_proportions_zip    <- file.path(local, "qmp_genera_abundances_BMIS.tsv.zip")
    scale_zip                   <- file.path(local, "41586_2020_2269_MOESM3_ESM.xlsx.zip")
    metadata_zip                <- file.path(local, "SraRunTable (34).csv.zip")
    motus_zip                   <- file.path(local, "PRJEB37249_motus_merged.tsv.zip")
    metaphlan4_zip              <- file.path(local, "PRJEB37249_MetaPhlAn_merged.tsv.zip")

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

    # ----- original proportions --

    if (file.exists(original_proportions_zip)) {
        zinfo <- unzip(original_proportions_zip, list = TRUE)
        tsv_name <- zinfo$Name[grepl("\\.tsv$", zinfo$Name)][1]
        if (is.na(tsv_name)) stop("No .tsv file found inside ", original_proportions_zip)
        tmpdir <- tempfile("orig_prop_")
        dir.create(tmpdir)
        unzip(original_proportions_zip, files = tsv_name, exdir = tmpdir)
        tsv_path <- file.path(tmpdir, tsv_name)
        counts_original <- read.delim(tsv_path, row.names = 1, check.names = FALSE)
        counts_original <- as.data.frame(counts_original)
        counts_original <- counts_original %>%
            rename(Sample = SampleID)
        counts_original <- counts_original %>%
            mutate(
            Sample = if_else(
                str_detect(Sample, "^x"),
                paste0("M0", Sample),
                Sample
            )
        )
        proportions_original <- sweep(counts_original, 2, colSums(counts_original), FUN = "/")
        old_colnames <- colnames(proportions_original)
        taxon_ids <- paste0("Taxon_", seq_along(old_colnames))
        colnames(proportions_original) <- taxon_ids
        tax_original <- data.frame(
        Taxon = taxon_ids,
        Genera = old_colnames,
        stringsAsFactors = FALSE
        )
        unlink(tmpdir, recursive = TRUE)
    }

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


    # ----- Metadata -----
    metadata <- NA
    if (file.exists(metadata_zip)) {
        metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
        metadata_con <- unz(metadata_zip, metadata_csv)
        metadata <- read.csv(metadata_con, row.names = "Submitter_id") %>%
        as.data.frame() %>%
        rownames_to_column("Sample")
    }

    # ----- Scale ----- 
    scale <- NA
    if (file.exists(scale_zip)) {
        zinfo     <- unzip(scale_zip, list = TRUE)
        xlsx_name <- zinfo$Name[grepl("\\.xlsx$", zinfo$Name)][1]
        if (is.na(xlsx_name)) stop("No .xlsx file in ", scale_zip)
        tmpdir     <- tempfile("scale_xlsx_")
        dir.create(tmpdir)
        unzip(scale_zip, files = xlsx_name, exdir = tmpdir)
        scale_path <- file.path(tmpdir, xlsx_name)
        COLUMN2 <- "Faecal microbial load (cells/g)"
        sheet2_full <- read_excel(scale_path, sheet = 3) %>% 
            as.data.frame()
        sheet2_full <- sheet2_full %>%
            rename(Sample = SampleID)
        sheet2_full <- sheet2_full %>%
            mutate(
            Sample = if_else(
                str_detect(Sample, "^x"),
                paste0("M0", Sample),
                Sample
            )
            )
        scale     <- sheet2_full %>% select(Sample, all_of(COLUMN2))
        metadata2 <- sheet2_full %>% select(-all_of(COLUMN2))
        unlink(tmpdir, recursive = TRUE)
    }

    metadata <- metadata %>%
        full_join(metadata2, by = "Sample")

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