require(stringr)
library(readxl)
library(dplyr)
library(R.utils)

parse_2023_pereira_nature_nervous <- function(rootPath = getwd(), paths = NULL) {
  localPath <- paste0(rootPath, "/2023_pereira_nature_nervous/")
  
  decompressed_file <- gunzip(paste0(localPath, "all-relevant-data.xlsx.gz"), remove = FALSE)
  
  totalCounts <- readxl::read_xlsx(decompressed_file, sheet = "Sup. Table 1", skip = 1)
  totalCounts <- totalCounts[, c("time (hours)", "condition", "medium", "Cells/mL")]
  totalCounts$index <- paste(totalCounts[["time (hours)"]], totalCounts$condition, totalCounts$medium, sep = "_")
  totalCounts$replicate <- ave(totalCounts$index, totalCounts$index, FUN = seq_along)
  totalCounts$index <- paste(totalCounts$index, totalCounts$replicate, sep = "_")
  totalCounts <- totalCounts[, c("index", "Cells/mL")]
  
  sampleMetadata <- readxl::read_xlsx(decompressed_file, sheet = "Sup. Table 2", skip = 3)
  
  # Step 1: Process sampleMetadata
  sampleMetadata <- sampleMetadata %>%
    mutate(index = paste(`time (hours)`, condition, medium, sep = "_")) %>%
    group_by(index) %>%
    mutate(replicate = seq_along(index)) %>%
    ungroup() %>%
    mutate(index = paste(index, replicate, sep = "_"))
  
  # Step 2: Intersect with totalCounts
  scale_data <- sampleMetadata %>%
    inner_join(totalCounts, by = "index") %>%
    dplyr::select("Sequencing sample ID", `Cells/mL`)
  
  # Step 3: Create scale matrix
  scale <- scale_data %>%
    pivot_wider(names_from = "Sequencing sample ID", values_from = `Cells/mL`) %>%
    as.matrix()
  
  relativeAbundances <- readxl::read_xlsx(decompressed_file, sheet = "Sup. Table 3", skip = 3)
  absoluteAbundances <- readxl::read_xlsx(decompressed_file, sheet = "Sup. Table 4", skip = 3)
  
  absoluteAbundances <- absoluteAbundances %>%
    rename(ASV_ID = `...1`)
  
  relativeAbundances <- relativeAbundances %>%
    rename(ASV_ID = `...1`)
  
  counts <- absoluteAbundances %>%
    dplyr::select(-c(Phylum, Class, Order, Family, Genus)) %>%
    column_to_rownames("ASV_ID") %>%
    as.matrix()
  
  proportions <- relativeAbundances %>%
    dplyr::select(-c(Phylum, Class, Order, Family, Genus)) %>%
    column_to_rownames("ASV_ID") %>%
    as.matrix()
  
  metadata <- sampleMetadata %>%
    column_to_rownames("Sequencing sample ID")
  
  tax <- absoluteAbundances %>%
    dplyr::select(ASV_ID, Phylum, Class, Order, Family, Genus) %>%
    column_to_rownames("ASV_ID")
  
  unlink(decompressed_file)
  
  return(list(counts = counts, scale = scale, proportions = proportions, tax = tax, metadata = metadata))
}