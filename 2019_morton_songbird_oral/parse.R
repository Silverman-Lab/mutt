library(tidyverse)

parse_2019_morton_songbird_oral <- function() {
  ## Read Counts
  counts <- readRDS("2019_morton_songbird_oral_counts.RDS")
  proportions <- apply(counts, 2, function(col) col/sum(col))

  ## Taxonomy Information
  raw_tax <- read.csv("taxonomy.tsv", sep="\t")
  tax <- raw_tax %>%
    mutate(Taxon=str_trim(Taxon)) %>%
    separate(
      Taxon,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = "\\s*;\\s*",
      extra = "drop",
      fill = "right"
  )
  row.names(tax) <- tax$Feature.ID
  tax <- tax[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

  ## Get metadata, add average flowcount, round to avoid fraction of a cell
  metadata <- data.frame(read.csv("oral_trimmed_metadata.csv", sep="\t", row.names=1))
  metadata$avg_norm_flowcounts <- round(rowMeans(metadata[,c("flow.cells.ul.1", "flow.cells.ul.2")]))
  scale <- metadata$avg_norm_flowcounts
  names(scale) <- row.names(metadata)

  return(list(counts=counts,
              proportions=proportions,
              scale=scale,
              metadata=metadata,
              tax=tax))
}

