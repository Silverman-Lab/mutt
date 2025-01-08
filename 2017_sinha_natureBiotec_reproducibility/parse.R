library(readxl)
parse_2024_sinha_natureBiotec_reproducibility <- function(paths = NULL) {
  local <- file.path("2017_sinha_natureBiotec_reproducibility/mbqc_org_files")
  
  dat <- read_tsv(paste0(local, "/mbqc_integrated_otus.tsv"))
  metadata <- dat[1:71,]
  metadata <- t(metadata)
  colnames(metadata) <- metadata[1,]
  metadata <- as.data.frame(metadata[-1,])
  
  # dat has metadata (row1-71) and otu counts (row72~) stitched together
  counts <- t(dat[72:nrow(dat),])
  colnames(counts) <- counts[1,]
  counts <- counts[-1,] # remove column names from row 1
  
  scale <- as.matrix(metadata$nt_count, ncol = 1)
  rownames(scale) <- rownames(metadata)
  
  # Extract taxonomy and the numerical IDs from count table
  tax_list <- strsplit(colnames(counts), "\\|")
  tax <- do.call(rbind, lapply(tax_list, function(taxa) {
    # Extract each level of taxonomy
    kingdom <- sub("k__", "", taxa[1])
    phylum <- sub("p__", "", taxa[2])
    class <- sub("c__", "", taxa[3])
    order <- sub("o__", "", taxa[4])
    family <- sub("f__", "", taxa[5])
    genus <- sub("g__", "", taxa[6])
    species <- sub("s__", "", taxa[7])
    
    # Extract the numerical ID (last part of each column names)
    id <- sub("^.*\\|", "", taxa[length(taxa)])
    
    c(kingdom, phylum, class, order, family, genus, species, id)
  }))
  
  # Convert to data frame and set column names
  tax <- as.data.frame(tax, stringsAsFactors = FALSE)
  colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ID")
  
  # Use the numerical IDs as rownames
  rownames(tax) <- tax$ID
  tax <- tax[, -ncol(tax)]  # Remove the "ID" column
  
  colnames(counts) <- rownames(tax)
  
  return(list(scale=scale, metadata=metadata, counts=t(counts), tax=tax))
}
