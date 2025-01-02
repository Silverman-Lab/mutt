library(readxl)
library(stringr)
parse_2020_zemb_microbiologyOpen_spike <- function(paths = NULL) {
  local <- file.path("~/Research/datarepo","2020_zemb_microbiologyOpen_spike/")
  local <- file.path("/storage/work/wpg5129/sparse_project/2020_zemb_microOpen_spike/")
  setwd(local)
  
  metadata <- read.csv(paste0(local,"zemb_metadata.csv"), row.names = 1)
  
  counts <- read.csv(paste0(local,"zemb_counts.csv"), row.names = 1)
  
  scale <- read.csv(paste0(local,"zemb_qPCR.csv"), row.names = 1)
  
  colnames(counts) <- paste0("oz1802-", unlist(lapply(str_split(colnames(counts), "Tube"), function(x) {x[2]})))
  counts <- counts[,colnames(counts)%in%rownames(metadata)]
  metadata <- metadata[rownames(metadata)%in%colnames(counts),]
  scale <- scale[rownames(scale)%in%colnames(counts),]
  
  all.equal(rownames(metadata), rownames((scale)))
  all.equal(colnames(counts), rownames(scale))
  
  return(list(scale=scale, metadata=metadata, counts=counts))
}
