require(stringr)
library(readxl)
parse_2024_jin_natureComm_technicalReplicates <- function(paths = NULL) {
  local <- "2022_jin_natureComm_technicalReplicates/"
  
  dat <- readxl::read_xlsx(paste0(local, "Supplementary Data 8.xlsx"), sheet = "Sequencing-determined counts", skip = 1)
  counts <- as.matrix(dat[,2:56])
  rownames(counts) <- dat$`cOTU-ID`
  tax <- dat[,57:62]
  rownames(tax) <- dat$`cOTU-ID`
  scale <- readxl::read_xlsx(paste0(local, "Supplementary Data 8.xlsx"), sheet = "Absolute total abundance", skip = 1)
  abs_counts <- readxl::read_xlsx(paste0(local, "Supplementary Data 8.xlsx"), sheet = "Absolute abundance", skip = 1)
  
  #extract metadata from sample id
  metadata <- as.data.frame(matrix(NA, nrow = nrow(scale), ncol = 5))
  rownames(metadata) <- scale$Sample
  colnames(metadata) <- c("diet", "subject", "location", "sample_type", "technical_replicate")
  
  #first two letter indicate diet type (VD, VS, CE)
  metadata$diet <- substr(scale$Sample, start = 1, stop = 2)
  
  #third letter indicate subject id within group. Note; subject in different group are different mice (i.e., VDa != VSa)
  metadata$subject <- substr(scale$Sample, start = 3, stop = 3)
  
  #sample collected in two different location (distal, proximal)
  metadata$location <- substr(scale$Sample, start = 4, stop = 7)
  
  #sample measured using either filter-residue (cell-sample) and flow-through (ecDNA-sample)
  sampleid <- rownames(metadata)
  metadata$sample_type <- unlist(lapply(strsplit(sampleid, split = "-"), function(x) x[2]))
  metadata$sample_type[is.na(metadata$sample_type)] <- "cell"
  
  #technical replicate if applicable (not every diet group have technical replicates)
  metadata$technical_replicate <- as.numeric(substr(scale$Sample, start = 8, stop = 8))
  
  #remove samples that doesn't match
  counts <- counts[,colnames(counts) %in% rownames(metadata)]
  metadata <- metadata[rownames(metadata) %in% colnames(counts),]
  scale <- scale[scale$Sample %in% colnames(counts),]
  abs_counts <- abs_counts[,colnames(abs_counts) %in% colnames(counts)]
  
  all.equal(colnames(counts), rownames(metadata))
  
  
  return(list(scale=scale, metadata=metadata, counts=counts, abs_counts = abs_counts, tax=tax))
}


