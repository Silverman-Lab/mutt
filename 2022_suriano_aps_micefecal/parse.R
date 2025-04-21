parse_2022_suriano_aps_micefecal <- function() {
  # Check for required packages
  required_pkgs <- c("tidyverse", "readxl", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
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
  
  metadata_zip <- "2022_suriano_aps_micefecal/SraRunTable.csv.zip"
  metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]  # list file inside zip
  metadata_path <- unzip(metadata_zip, files = metadata_csv, exdir = tempdir(), overwrite = TRUE)
  metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

  metadata$Sample_name <- gsub("^(.*?_D)(\\d{2})\\.(\\d+)$", "D\\2-\\3", metadata$Sample_name)
  rownames(metadata) <- metadata$Sample_name
  metadata <- subset(metadata, select = -Sample_name)

  ## Read scale sheet
  scale_zip <- "2022_suriano_aps_micefecal/Supplementary Table 6.xlsx.zip"
  scale_xlsx <- unzip(scale_zip, list = TRUE)$Name[1]
  scale_path <- unzip(scale_zip, files = scale_xlsx, exdir = tempdir(), overwrite = TRUE)
  
  scale <- read_excel(scale_path, sheet = "Microbial loads")
  scale <- scale %>%
    select(Sample, `Cells.g.of.fecal.sample`) %>%
    pivot_wider(names_from = Sample, values_from = `Cells.g.of.fecal.sample`)
  
  return(list(
    # counts = counts,
    # proportions = proportions,
    # tax = tax,
    scale = scale,
    metadata = metadata
  ))
}