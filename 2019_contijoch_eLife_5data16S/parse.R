library(readxl)
library(biomformat)
library(readr)
parse_2019_contijoch_eLife_5data16S <- function(paths = NULL) {
  local <- file.path("2019_contijoch_eLife_5data16S/")
  
  dat <- read_biom(paste0(local,"16S_full_combined_OTU_table.biom"))
  
  counts <- as(biom_data(dat), "matrix")
  
  # data associated with IBD (Crohn's disease) : 16S rRNA sequencing data (data 1-3)
  metadata1.1 <- read_xlsx(paste0(local,"elife-40553-fig1-data1-v1.xlsx"), skip = 2)
  metadata1.1$SampleID <- str_replace_all(metadata1.1$SampleID, "-|_", ".")
  dim(metadata1.1)
  sum(metadata1.1$SampleID %in% colnames(counts)) # one missing (108/109)
  # drop samples don't exist in metadata
  counts1.1 <- counts[,colnames(counts) %in% metadata1.1$SampleID]
  # drop samples don't exist in counts
  metadata1.1 <- metadata1.1[metadata1.1$SampleID %in% colnames(counts1.1),]
  # reorder counts columns (same order as metadata)
  counts1.1 <- counts1.1[,metadata1.1$SampleID]
  
  # 1 DNA molecule weighs 10^-15 grams
  # microbiota density is defined as the total DNA extracted from each sample (in mg) per mg of freshsample.
  scale1.1 <- data.frame(total = metadata1.1$`Sample Mass` * (metadata1.1$`Microbiota Density` * 10^-6) / (10^-15))
  rownames(scale1.1) <- metadata1.1$SampleID
  
  metadata1.3 <- read_xlsx(paste0(local,"elife-40553-fig1-data3-v1.xlsx"), skip = 2)
  metadata1.3$SampleID <- str_replace_all(metadata1.3$SampleID, "-|_", ".")
  dim(metadata1.3)
  sum(metadata1.3$SampleID %in% colnames(counts)) # all match
  # drop samples don't exist in metadata
  counts1.3 <- counts[,colnames(counts) %in% metadata1.3$SampleID]
  # drop samples don't exist in counts
  metadata1.3 <- metadata1.3[metadata1.3$SampleID %in% colnames(counts1.3),]
  # reorder counts columns (same order as metadata)
  counts1.3 <- counts1.3[,metadata1.3$SampleID]
  scale1.3 <- data.frame(total = metadata1.3$`Sample Mass` * (metadata1.3$`Microbiota Density` * 10^-6) / (10^-15))
  rownames(scale1.3) <- metadata1.3$SampleID
  
  metadata2.1 <- read_xlsx(paste0(local,"elife-40553-fig2-data1-v1.xlsx"), skip = 2)
  metadata2.1$SampleID <- str_replace_all(metadata2.1$SampleID, "-|_", ".")
  dim(metadata2.1)
  sum(metadata2.1$SampleID %in% colnames(counts)) # 385/649 match
  # drop samples don't exist in metadata
  counts2.1 <- counts[,colnames(counts) %in% metadata2.1$SampleID]
  # drop samples don't exist in counts
  metadata2.1 <- metadata2.1[metadata2.1$SampleID %in% colnames(counts2.1),]
  # reorder counts columns (same order as metadata)
  counts2.1 <- counts2.1[,metadata2.1$SampleID]
  scale2.1 <- data.frame(total = metadata2.1$`Sample Mass` * (metadata2.1$`Microbiota Density` * 10^-6) / (10^-15))
  rownames(scale2.1) <- metadata2.1$SampleID
  
  metadata4.1 <- read_xlsx(paste0(local,"elife-40553-fig4-data1-v1.xlsx"), skip = 2)
  metadata4.1$SampleID <- str_replace_all(metadata4.1$SampleID, "-|_", ".")
  dim(metadata4.1)
  sum(metadata4.1$SampleID %in% colnames(counts)) # all match
  # drop samples don't exist in metadata
  counts4.1 <- counts[,colnames(counts) %in% metadata4.1$SampleID]
  # drop samples don't exist in counts
  metadata4.1 <- metadata4.1[metadata4.1$SampleID %in% colnames(counts4.1),]
  # reorder counts columns (same order as metadata)
  counts4.1 <- counts4.1[,metadata4.1$SampleID]
  scale4.1 <- data.frame(total = metadata4.1$`Sample Mass` * (metadata4.1$`Microbiota Density` * 10^-6) / (10^-15))
  rownames(scale4.1) <- metadata4.1$SampleID
  
  metadata4.2 <- read_xlsx(paste0(local,"elife-40553-fig4-data2-v1.xlsx"), skip = 2)
  metadata4.2$SampleID <- str_replace_all(metadata4.2$SampleID, "-|_", ".")
  dim(metadata4.2)
  sum(metadata4.2$SampleID %in% colnames(counts)) # 4 missing
  # drop samples don't exist in metadata
  counts4.2 <- counts[,colnames(counts) %in% metadata4.2$SampleID]
  # drop samples don't exist in counts
  metadata4.2 <- metadata4.2[metadata4.2$SampleID %in% colnames(counts4.2),]
  # reorder counts columns (same order as metadata)
  counts4.2 <- counts4.2[,metadata4.2$SampleID]
  scale4.2 <- data.frame(total = metadata4.2$`Sample Mass` * (metadata4.2$`Microbiota Density` * 10^-6) / (10^-15))
  rownames(scale4.2) <- metadata4.2$SampleID
  
  # tax - note from author ; 
  # OTU ID in these samples corresponds to those from the (now no longer maintained) greengenes database gg_13_8, clustered at 97% similarity,
  # note; gg_13_8_99 is used instead (couldn't find gg_13_8_97)
  
  # Load the gg.tax file
  gg_tax <- read.delim(paste0(local, "gg_13_8_99.gg.tax"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Add column names (adjust as needed based on file structure)
  colnames(gg_tax) <- c("OTU_ID", "Taxonomy")
  head(gg_tax)
  
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  # Split strings into components
  gg_split <- strsplit(gg_tax$Taxonomy, "\\;")
  
  # Create a matrix with 7 columns
  gg_matrix <- do.call(rbind, lapply(gg_split, function(entry) {
    # Extract ranks, ensuring a length of 7 (filling with NA if missing)
    ranks <- sapply(ranks, function(rank) {
      matched <- grep(paste0("^", substr(rank, 1, 1), "__"), entry, value = TRUE)
      if (length(matched) == 0) return(NA) else return(sub("^[a-z]__", "", matched))
    })
    return(ranks)
  }))
  
  # Convert to data frame for better visualization
  tax <- as.data.frame(gg_matrix, stringsAsFactors = FALSE)
  colnames(tax) <- ranks
  
  #all.equal(colnames(counts4.2), metadata4.2$SampleID)
  #all.equal(colnames(counts4.2), rownames(scale4.2))
  
  out <- list()
  out[[1]] <- list(scale=scale1.1, metadata=metadata1.1, counts=counts1.1, tax=tax)
  out[[2]] <- list(scale=scale1.3, metadata=metadata1.3, counts=counts1.3, tax=tax)
  out[[3]] <- list(scale=scale2.1, metadata=metadata2.1, counts=counts2.1, tax=tax)
  out[[4]] <- list(scale=scale4.1, metadata=metadata4.1, counts=counts4.1, tax=tax)
  out[[5]] <- list(scale=scale4.2, metadata=metadata4.2, counts=counts4.2, tax=tax)
  names(out) <- paste0("data", c(1.1, 1.3, 2.1, 4.1, 4.2))
  
  return(list(list(scale=scale, metadata=metadata, counts=counts, tax=tax)))
}
