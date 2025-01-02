library(stringr)
library(readxl)
parse_2024_jin_pnas_semen <- function(paths = NULL) {
  local <- file.path("/storage/work/wpg5129/sparse_project","2024_jin_pnas_semen")
  
  metadata <- as.data.frame(read_table(paste0(local, "/metadata4.txt")))
  rownames(metadata) <- metadata$`sample-ID`
  
  counts <- as.data.frame(read_tsv(paste0(local, "/table.tsv"), skip = 1))
  rownames(counts) <- counts$`#OTU ID`
  
  tax <- read_tsv(paste0(local, "/taxonomy.tsv"))
  rownames(tax) <- tax$`Feature ID`
  
  scale <- read_table(paste0(local, "/CFU_data.txt"))
  rownames(scale) <- scale$`sample-ID`
  
  ## TODO some samples do not match between metadata and counts
  counts <- counts[,colnames(counts) %in% rownames(metadata)]
  metadata <- metadata[rownames(metadata) %in% colnames(counts),]
  
  all.equal(colnames(counts), rownames(metadata))
  
  return(list(scale=scale, metadata=metadata, counts=counts, tax=tax))
}

