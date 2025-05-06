parse_2022_zaramela_msystems_synDNA <- function(raw = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2022_zaramela_msystems_synDNA")

    # ----- File paths -----
    motus_zip      <- file.path(local, "PRJNA940499_motus_merged.tsv.zip")
    metaphlan4_zip <- file.path(local, "PRJNA940499_MetaPhlAn_merged.tsv.zip")
    metadata_SRA_zip     <- file.path(local, "SraRunTable (34).csv.zip")
    metadatareadin       <- file.path(local, "synDNA_metadata_updated.tsv.zip")
    countsreadin         <- file.path(local, "HMP_frequency_table.tsv.zip")
    scalereadin          <- file.path(local, "Saliva_Flow_data.tsv.zip")
    taxreadin            <- file.path(local, "taxaID_length_fulllineage.tsv.zip")

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

    # ----- scale and metadata -----

    metadata1 <- read_zipped_table(metadata_SRA_zip, row.names=NULL) %>%
                    rename(Accession = Run, Sample = `Sample Name`)
    metadata2 <- read_zipped_table(metadatareadin, sep="\t")
    metadata2 <- metadata2 %>% rename(SampleID = Sample) %>% mutate(Sample = paste(SampleID, Pool, sep = "_"))
    metadata = full_join(metadata1, metadata2, by="Sample")
    scale <- read_zipped_table(scalereadin, row.names=NULL, sep="\t") %>%
             rename(Sample = SampleID)
    metadata = left_join(metadata, scale %>% select(c("Sample", "gender", "saliva_weight_g", 
                "saliva_volume_mL_collected_in_5_min", "saliva_flow_rate_mL_per_min")), by="Sample")
    scale = scale %>% select(-c("gender")) 

    # ----- counts and proportions and taxa -----

    counts <- read_zipped_table(countsreadin, row.names=NULL, sep="\t") %>%
                group_by(OTUID) %>%  
                summarise(across(where(is.numeric), sum, na.rm = TRUE))
    counts <- as.data.frame(counts)
    rownames(counts) <- counts$OTUID
    counts$OTUID <- NULL
    counts = as.data.frame(t(counts))

    tmp_file <- unzip(taxreadin, exdir = tempdir())[1]
    tax <- read_tsv(tmp_file, quote = "", col_types = cols(.default = "c")) %>%
    rename(Kingdom = superkingdom, Phylum = phylum, Class = class, 
            Order = order, Family = family, Genus = genus, Species = species, 
            OTUID = GenomeID)

    tax = make_taxa_label(tax)
    rownames(tax) = tax$OTUID

    if (!raw) {
        matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
        colnames(counts) <- matched_taxa
        counts <- as.data.frame(t(rowsum(t(counts), group = colnames(counts))))
    }

    # --- Compute proportions from counts ---
    proportions <- sweep(counts, MARGIN = 1,STATS  = rowSums(counts), FUN = "/")

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
