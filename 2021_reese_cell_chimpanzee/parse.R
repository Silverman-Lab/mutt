parse_2021_reese_cell_chimpanzee <- function() {
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

    metadata_zip <- "2021_reese_cell_chimpanzee/SraRunTable.csv.zip"
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1] # list file inside zip
    metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
    metadata <- as.data.frame(read.csv(metadata_path, stringsAsFactors = FALSE))

    rownames(metadata) <- metadata$Sample_name
    metadata <- subset(metadata, select = -Sample_name)
    rownames(metadata) <- gsub("c$", "", rownames(metadata))

    ## Read scale sheet
    scale_zip <- "2021_reese_cell_chimpanzee/1-s2.0-S0960982220316523-mmc3.xlsx.zip"
    scale_xlsx <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_xlsx, exdir = tempdir(), overwrite = TRUE)

    scale <- as.data.frame(read_excel(scale_path, sheet = "S2 Metadata"))
    rownames(scale) <- scale$`Sample ID`
    scale <- scale[rownames(metadata), , drop = FALSE]


    scale_data <- scale %>%
        select(`Sample ID`, `16S copies per g feces`) %>%
        pivot_wider(names_from = `Sample ID`, values_from = `16S copies per g feces`)

    scale_extra <- scale %>%
        select(-`Sample ID`, -`16S copies per g feces`)

    scale_extra <- scale_extra[rownames(metadata), , drop = FALSE]

    metadata <- cbind(metadata, scale_extra)

    return(list(
        # counts = counts,
        # proportions = proportions,
        # tax = tax,
        scale = scale,
        metadata = metadata
    ))
}
