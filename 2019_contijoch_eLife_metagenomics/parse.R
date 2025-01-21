library(biomformat)
library(readr)
parse_2019_contijoch_eLife_metagenomic <- function(paths = NULL) {
  local <- file.path("~/2019_contijoch_eLife_5data16S/")
  
  # data associated with FMT in patients with rCDI : metagenomic data (data 4) n=15
  metadata <- read_table(paste0(local,"metagenomics_sample_metadata.txt"))

  counts <- as.data.frame(read_table(paste0(local, "metagenomics_combined_full.txt")))
  
  # Define the ranks in order
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")

  # Split strings into components
  tax_split <- strsplit(counts$Taxa, "\\|")

  # Create a matrix with 7 columns
  tax_matrix <- do.call(rbind, lapply(tax_split, function(entry) {
    # Extract ranks, ensuring a length of 7 (filling with NA if missing)
    ranks <- sapply(ranks, function(rank) {
      matched <- grep(paste0("^", substr(rank, 1, 1), "__"), entry, value = TRUE)
      if (length(matched) == 0) return(NA) else return(sub("^[a-z]__", "", matched))
    })
    return(ranks)
  }))

  # Convert to data frame for better visualization
  tax <- as.data.frame(tax_matrix, stringsAsFactors = FALSE)
  colnames(tax) <- ranks
  rownames(tax) <- paste0("seq", 1:nrow(counts))
  
  rownames(counts) <- rownames(tax)
  counts <- counts[,-1]
  
  scale <- data.frame(total = metadata$sample_mass_mg * (metadata$sample_microbial_density_ug_per_mg * 10^-6) / (10^-15))
  rownames(scale) <- metadata$Sample_ID
  
  all.equal(colnames(counts), metadata$Sample_ID)
  all.equal(rownames(scale), metadata$Sample_ID)

  return(list(scale=scale, metadata=metadata, counts=counts, tax=tax))
}
