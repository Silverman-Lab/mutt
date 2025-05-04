parse_2021_reese_cell_chimpanzee <- function(raw = FALSE) {
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

    # ---------- local paths --------------------
    local               <- file.path("2021_reese_cell_chimpanzee")

    # ---------- file paths ---------------------
    repro_counts_rds_zip<- file.path(local, "PRJEB39807_dada2_counts.rds.zip")
    repro_tax_zip       <- file.path(local, "PRJEB39807_dada2_taxa.rds.zip")
    metadata_zip        <- file.path(local, "SraRunTable.csv.zip")
    scale_zip           <- file.path(local,"1-s2.0-S0960982220316523-mmc3.xlsx.zip")

    # ---- scale and metadata -------------------
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1] # list file inside zip
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- as.data.frame(read.csv(metadata_path, stringsAsFactors = FALSE))
    rownames(metadata) <- metadata$Sample_name
    metadata <- subset(metadata, select = -Sample_name)
    rownames(metadata) <- gsub("c$", "", rownames(metadata))
    scale_xlsx <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_xlsx, exdir = tempdir(), overwrite = TRUE)
    scale <- as.data.frame(read_excel(scale_path, sheet = "S2 Metadata"))
    rownames(scale) <- scale$`Sample ID`
    scale <- scale[rownames(metadata), , drop = FALSE]
    scale_data <- scale %>%
        select(`Sample ID`, `16S copies per g feces`) %>%
        rename(Sample_name = `Sample ID`)

    scale_extra <- scale %>% select(-`Sample ID`, -`16S copies per g feces`)
    scale_extra <- scale_extra[rownames(metadata), , drop = FALSE]

    metadata <- cbind(metadata, scale_extra) %>%
        rename(
            Sample_name = Submitter_Id,
            Accession = Run
        )

    scale <- scale_data %>%
        left_join(metadata %>% select(Accession, Sample_name), by = "Sample_name")


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
