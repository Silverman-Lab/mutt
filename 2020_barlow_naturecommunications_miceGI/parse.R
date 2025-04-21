parse_2020_zemb_microOpen_spike <- function() {
    # Check for required packages
    required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
            "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function."
        )
    }

    # Load libraries
    library(tidyverse)
    library(readxl)
    library(stringr)
    library(readr)

    # --- Metadata ---
    metadata_zip <- "2020_barlow_naturecommunications_miceGI/metadata.xlsx.zip"
    metadata_xlsx <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_path <- unzip(metadata_zip, files = metadata_xlsx, exdir = tempdir(), overwrite = TRUE)
    metadata <- read_excel(metadata_path, sheet="Sheet1")

    sra_metadata_zip <- "2020_zemb_microOpen_spike/SraRunTable.csv.zip"
    sra_metadata_csv <- unzip(sra_metadata_zip, list = TRUE)$Name[1]
    sra_metadata_path <- unzip(sra_metadata_zip, files = sra_metadata_csv, exdir = tempdir(), overwrite = TRUE)
    sra_metadata <- read.csv(sra_metadata_path, row.names = 1, stringsAsFactors = FALSE)
    rownames(sra_metadata) <- sra_metadata$Sample.Name
    sra_metadata <- subset(sra_metadata, select = -Sample.Name)
    rownames(sra_metadata) <- paste0("oz1802-", str_extract(rownames(sra_metadata), "(?<=Tube).*")) 

    # --- Counts ---
    counts_zip <- "2020_zemb_microOpen_spike/zemb_counts.csv.zip"
    counts_csv <- unzip(counts_zip, list = TRUE)$Name[1]
    counts_path <- unzip(counts_zip, files = counts_csv, exdir = tempdir(), overwrite = TRUE)
    counts <- read.csv(counts_path, row.names = 1, stringsAsFactors = FALSE)

    # Rename count columns
    colnames(counts) <- paste0("oz1802-", str_extract(colnames(counts), "(?<=Tube).*"))

    # --- Scale (qPCR) ---
    scale_zip <- "2020_zemb_microOpen_spike/zemb_qPCR.csv.zip"
    scale_csv <- unzip(scale_zip, list = TRUE)$Name[1]
    scale_path <- unzip(scale_zip, files = scale_csv, exdir = tempdir(), overwrite = TRUE)
    scale <- read.csv(scale_path, row.names = 1, stringsAsFactors = FALSE)

    counts <- counts[, colnames(counts) %in% rownames(metadata)]
    metadata <- metadata[rownames(metadata) %in% colnames(counts), ]
    scale <- scale[rownames(scale) %in% colnames(counts), ]
    sra_metadata <- sra_metadata[rownames(sra_metadata) %in% colnames(counts), ]

    # Merge sra_metadata and metadata
    metadata <- cbind(metadata, sra_metadata[rownames(metadata), , drop = FALSE])

    # Check alignment
    all.equal(rownames(metadata), rownames(scale))
    all.equal(colnames(counts), rownames(scale))

    # Return structured list
    return(list(
        counts = counts,
        metadata = metadata,
        scale = scale
    ))
}
