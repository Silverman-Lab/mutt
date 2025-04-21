parse_2023_kallastu_research_foodscience_food <- function() {
    # Check for required packages
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

    # ## Read Counts
    # counts <- NA
    # raw_tax <- row.names(counts)
    # row.names(counts) <- paste0("Taxon_", 1:nrow(counts))
    # proportions <- apply(counts, 2, function(col) col/sum(col))
    #
    # ## Create Taxa
    # raw_tax <- data.frame(taxa=raw_tax)
    # tax <- raw_tax %>%
    #   mutate(taxa=str_trim(taxa)) %>%
    #   separate(
    #     taxa,
    #     into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    #     sep = "\\s*;\\s*",
    #     extra = "drop",
    #     fill = "right"
    # )
    # row.names(tax) <- row.names(counts)

    metadata_zip <- "2023_kallastu_research_foodscience_food/SraRunTable.csv.zip"
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1] # list file inside zip
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- as.data.frame(read.csv(metadata_path, stringsAsFactors = FALSE))

    rownames(metadata) <- metadata$`Sample.Name`
    metadata <- subset(metadata, select = -`Sample.Name`)

    ## Scale is directly presented as a table in the paper
    scale_data <- data.frame(
        "20St" = c(1.55e9, NA, 1.56e9, 1.16e9),
        `K6 PMA` = c(1.28e8, 1.53e7, 5.96e7, 1.71e8),
        `K6 TOT` = c(1.86e9, NA, 1.26e8, 3.06e8),
        `K5 PMA` = c(3.82e8, 3.64e7, 7.04e7, 3.73e8),
        `K5 TOT` = c(5.90e8, NA, 1.34e8, 5.66e8),
        `K4 PMA` = c(3.97e8, 1.17e8, 2.61e7, 4.12e8),
        `K4 TOT` = c(8.56e8, NA, 1.49e8, 4.60e8),
        `20St 2.5PMA` = c(1.51e9, NA, NA, NA),
        `20St 1PMA`   = c(1.39e9, NA, NA, NA),
        `20St 0.5PMA` = c(1.74e9, NA, NA, NA)
    )
    rownames(scale_data) <- c("Spike_In", "Plating_cfu_per_g", "FC", "qPCR")
    colnames(scale_data) <- c(
        "20St", "K6 PMA", "K6 TOT", "K5 PMA", "K5 TOT",
        "K4 PMA", "K4 TOT", "20St 2.5PMA", "20St 1PMA", "20St 0.5PMA"
    )

    return(list(
        # counts = counts,
        # proportions = proportions,
        # tax = tax,
        scale = scale,
        metadata = metadata
    ))
}