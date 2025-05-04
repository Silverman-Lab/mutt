parse_2020_regalado_isme_metagenomicsbacteriaandfungalsequencingqPCR <- function(raw = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2020_regalado_isme_metagenomicsbacteriaandfungalsequencingqPCR")

    # ----- File paths -----
    motus_zip            <- file.path(local, "PRJNA31530_motus_merged.tsv.zip")
    metaphlan4_zip       <- file.path(local, "PRJNA31530_MetaPhlAn_merged.tsv.zip")
    metadata_zip         <- file.path(local, "SraRunTable (38).csv.zip")
    ITS_metadata_zip     <- file.path(local, "ITS_blackblue_metadata_v2.txt.zip")
    scale_16s_zip         <- file.path(local, "")
    repro_counts_zips <- c(
    file.path(local, "PRJNA31530_dada2_counts.rds.zip")#,
    #file.path(local, "PRJNA31530_SILVA_counts.rds.zip")
    )

    repro_tax_zips <- c(
    file.path(local, "PRJNA31530_dada2_taxa.rds.zip")#,
    #file.path(local, "PRJNA31530_SILVA_taxa.rds.zip")
    )


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

    # ----- Initialize everything as NA -----
    counts_original_16s <- NA
    proportions_original_16s <- NA
    tax_original_16s <- NA
    counts_original_ITS <- NA
    proportions_original_ITS <- NA
    tax_original_ITS <- NA
    counts_original_meta <- NA
    proportions_original_meta <- NA
    tax_original_meta <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
    mOTU3_counts <- NA
    mOTU3_proportions <- NA
    mOTU3_tax <- NA
    MetaPhlAn4_counts <- NA
    MetaPhlAn4_proportions <- NA
    MetaPhlAn4_tax <- NA
    counts_ITS_reprocessed <- NA
    proportions_ITS_reprocessed <- NA
    tax_ITS_reprocessed <- NA

    # ------ original counts ------

    # ---- scale and metadata -----
    scale_16s_df     <- read_zipped_table(scale_16s_zip)
    scale_meta_df    <- read_zipped_table(scale_meta_zip)
    metadata_16s_df  <- read_zipped_table(metadata_16s_zip)
    metadata_meta_df <- read_zipped_table(metadata_meta_zip)

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

    return(list(
    counts = list(
                original = counts,
                reprocessed = list(
                    amplicon = counts_reprocessed,
                    shotgun = list(
                        mOTU3 = mOTU3_counts,
                        MetaPhlan4 = MetaPhlAn4_counts
                    ),
                    ITS = counts_ITS_reprocessed
                )

    ),
    proportions = list(
                original = proportions,
                reprocessed = list(
                    amplicon = proportions_reprocessed,
                    shotgun = list(
                        mOTU3 = mOTU3_proportions,
                        MetaPhlan4 = MetaPhlAn4_proportions
                    ),
                    ITS = proportions_ITS_reprocessed
                )
    ),
    tax = list(
                original = tax,
                reprocessed = list(
                    amplicon = tax_reprocessed,
                    shotgun = list(
                        mOTU3 = mOTU3_tax,
                        MetaPhlan4 = MetaPhlAn4_tax
                    ),
                    ITS = tax_ITS_reprocessed
                )

    ),
    scale = list(
                amplicon = scale_16s_df,
                shotgun = scale_meta_df,
                ITS = scale_ITS_df
    ),
    metadata = list(
                amplicon = metadata_16s_df,
                shotgun = metadata_meta_df,
                ITS = metadata_ITS_df
    )
    ))
}