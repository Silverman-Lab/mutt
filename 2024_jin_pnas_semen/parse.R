require(stringr)
library(readxl)
library(dplyr)

# removing the control samples as no scale available for them 
parse_2024_jin_pnas_semen <- function(rootPath = getwd(), paths = NULLL) {
  localPath <- paste0(rootPath, "/2024_jin_pnas_semen/")
  
  metadata <- as.data.frame(read_table(paste0(localPath, "metadata4.txt")))
  rownames(metadata) <- metadata$`sample-ID`
  metadata <- subset(metadata, select = -`sample-ID`)
  metadata <- metadata[!(rownames(metadata) %in% c("S1", "S13", "SN", "SNTC")), ]
  
  
  counts <- as.data.frame(read_tsv(paste0(localPath, "table.tsv"), skip = 1))
  rownames(counts) <- counts$`#OTU ID`
  counts <- subset(counts, select = -`#OTU ID`)
  counts <- t(counts)
  counts <- counts[!(rownames(counts) %in% c("S1", "SN", "SNTC")), ]
  
  
  tax <- as.data.frame(read_tsv(paste0(localPath, "taxonomy.tsv")), skip = 1)
  rownames(tax) <- tax$`Feature ID`
  tax <- subset(tax, select = -`Feature ID`)
  
  scale <- as.data.frame(read_table(paste0(localPath, "CFU_data.txt")), skip = 1)
  rownames(scale) <- scale$`sample-ID`
  scale <- subset(scale, select = -`sample-ID`)
  scale <- t(scale)
  
  return(list(scale=scale, metadata=metadata, counts=counts, tax=tax))
}