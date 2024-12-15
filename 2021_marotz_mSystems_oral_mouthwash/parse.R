library(tidyverse)

parse_2021_marotz_mSystems_oral_mouthwash <- function() {
  ## Read Counts
  counts <- readRDS("2021_marotz_mSystems_oral_mouthwash.RDS")
  raw_tax <- row.names(counts)
  row.names(counts) <- paste0("Taxon_", 1:nrow(counts))
  proportions <- apply(counts, 2, function(col) col/sum(col))

  ## Create Taxa Metadata
  raw_tax <- data.frame(taxa=raw_tax)
  tax <- raw_tax %>%
    mutate(taxa=str_trim(taxa)) %>%
    separate(
      taxa,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = "\\s*;\\s*",
      extra = "drop",
      fill = "right"
  )
  row.names(tax) <- row.names(counts)

  metadata <- data.frame(read.csv("T3_SRS_metadata_ms.txt", sep="\t", row.names=1))
  metadata <- metadata[colnames(counts),]
  scale <- metadata$FC_avg_cells_per_ul
  names(scale) <- row.names(metadata)

  return(list(counts=counts,
              proportions=proportions,
              scale=scale,
              metadata=metadata,
              tax=tax))
}

