parse_2023_pereira_nature_nervous <- function(raw = FALSE) {
    required_pkgs <- c("tidyverse", "readxl")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }
    
    require(stringr)
    library(tidyverse)

    # ---------- MANAN PROCESSED BELOW ----------

    # localPath <- file.path("/2023_pereira_nature_nervous/")

    # decompressed_file <- gunzip(paste0(localPath, "all-relevant-data.xlsx.gz"), remove = FALSE)
    
    # totalCounts <- readxl::read_xlsx(decompressed_file, sheet = "Sup. Table 1", skip = 1)
    # totalCounts <- totalCounts[, c("time (hours)", "condition", "medium", "Cells/mL")]
    # totalCounts$index <- paste(totalCounts[["time (hours)"]], totalCounts$condition, totalCounts$medium, sep = "_")
    # totalCounts$replicate <- ave(totalCounts$index, totalCounts$index, FUN = seq_along)
    # totalCounts$index <- paste(totalCounts$index, totalCounts$replicate, sep = "_")
    # totalCounts <- totalCounts[, c("index", "Cells/mL")]
    
    # sampleMetadata <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 2", skip = 3)
    
    # # Step 1: Process sampleMetadata
    # sampleMetadata <- sampleMetadata %>%
    #   mutate(index = paste(`time (hours)`, condition, medium, sep = "_")) %>%
    #   group_by(index) %>%
    #   mutate(replicate = seq_along(index)) %>%
    #   ungroup() %>%
    #   mutate(index = paste(index, replicate, sep = "_"))
    
    # # Step 2: Intersect with totalCounts
    # scale_data <- sampleMetadata %>%
    #   inner_join(totalCounts, by = "index") %>%
    #   dplyr::select("Sequencing sample ID", `Cells/mL`)
    
    # # Step 3: Create scale matrix
    # scale <- scale_data %>%
    #   pivot_wider(names_from = "Sequencing sample ID", values_from = `Cells/mL`) %>%
    #   as.matrix()
    
    # relativeAbundances <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 3", skip = 3)
    # absoluteAbundances <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 4", skip = 3)
    
    # absoluteAbundances <- absoluteAbundances %>%
    #   rename(ASV_ID = `...1`)
    
    # relativeAbundances <- relativeAbundances %>%
    #   rename(ASV_ID = `...1`)
    
    # countsDataFile <- unzip(paste0(localPath, "DADA2_counts_as_matrix_Drugs.txt.zip"))
    
    # countsData <- read.table(countsDataFile, header = TRUE, sep = "\t")
    
    # counts <- countsData %>% column_to_rownames("X") %>% as.matrix()
    
    # proportions <- relativeAbundances %>%
    #   dplyr::select(-c(Phylum, Class, Order, Family, Genus)) %>%
    #   column_to_rownames("ASV_ID") %>%
    #   as.matrix()
    
    # metadata <- sampleMetadata %>%
    #   column_to_rownames("Sequencing sample ID")
    
    # tax <- absoluteAbundances %>%
    #   dplyr::select(ASV_ID, Phylum, Class, Order, Family, Genus) %>%
    #   column_to_rownames("ASV_ID")
    
    # unlink(decompressed_file)
    # unlink(countsDataFile)

    # ------ MAXWELL PROCESSED BELOW ------------

    # ----- Local base directory -----
    local <- file.path("2023_pereira_nature_nervous")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJNA1033532_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA1033532_dada2_taxa.rds.zip")
    scale_16s_zip        <- file.path(local, "Pereira2023_scale.csv.zip")
    counts_16s_zip       <- file.path(local, "Pereira_2023_16S.csv.zip")
    metadata_16s_zip     <- file.path(local, "Pereira_2023_metadata.csv.zip")
    sra_zip              <- file.path(local, "SraRunTable (38).csv.zip")

    read_zipped_csv <- function(zip_path) {
        if (file.exists(zip_path)) {
            csv_file <- unzip(zip_path, list = TRUE)$Name[1]
            read.csv(unz(zip_path, csv_file), row.names = 1, check.names = FALSE)
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

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA

    # ------ original counts ------
    counts_original <- read_zipped_csv(counts_16s_zip) 

    if (!is.na(counts_original)[1]) {
        original_taxa <- colnames(counts_original)

        # Create taxa mapping data frame
        tax_original <- data.frame(
        Taxa = original_taxa,
        stringsAsFactors = FALSE
        )

        # ------ proportions from counts ------
        proportions_original <- sweep(counts_original, MARGIN = 1,STATS  = rowSums(counts_original), FUN = "/")

    } else {
        proportions_original <- NA
        tax_original <- NA
    }

    # ---- scale and metadata -----
    scale     <- read_zipped_csv(scale_16s_zip)
    metadata  <- read_zipped_csv(metadata_16s_zip)
    sra       <- read_zipped_csv(sra_zip) %>% rename(Accession = Run)

    # # ----- Reprocessed counts from RDS ZIP -----
    # temp_rds <- tempfile(fileext = ".rds")
    # unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

    # rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
    # if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
    # counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

    # # ----- Taxonomy reprocessed -----
    # temp_tax <- tempfile(fileext = ".rds")
    # unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

    # tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
    # if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
    # tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

    
    # # ----- Convert sequences to lowest rank taxonomy found and update key -----
    # tax_reprocessed = make_taxa_label(tax_reprocessed)

    # # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    # # accessions to sampleIDs is study specific: IF NEED BE

    # # taxa
    # if (!raw) {
    #     matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
    #     colnames(counts_reprocessed) <- matched_taxa
    #     counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
    # }

    # # proportions reprocessed
    # proportions_reprocessed = counts_reprocessed
    # proportions_reprocessed[-1] <- lapply(
    #     counts_reprocessed[-1],
    #     function(col) col / sum(col)
    # )

    if (!raw) {
        counts_original = fill_na_zero_numeric(counts_original)
        proportions_original = fill_na_zero_numeric(proportions_original)
        #counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        #proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }


    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original,
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = proportions_original,
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = tax_original,
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))
}