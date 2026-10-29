## Paper

This is from the songbird paper: https://www.nature.com/articles/s41467-019-10656-5

## Files

The raw files oral_trimmed_deblur.biom, oral_trimmed_metadata.csv, and taxonomy.tsv were downloaded from https://github.com/knightlab-analyses/reference-frames.

The file 2019_morton_songbird_oral_counts.RDS is the counts table in a better format than biom,
the code for going from biom to counts is:

```R
library(biomformat)
## Load sequence count data
biom_data <- read_hdf5_biom("oral_trimmed_deblur.biom")
Y <- do.call(rbind, biom_data$data)
row.names(Y) <- sapply(biom_data$rows, function(item){item$id})
```

## Notes on Data

* Metadata: the metadata contains many columns related to the scale including: flow.cell.5min.1, flow.cells.ul.1, and flow.cells.ul.2. The first two are not normalized by saliva sample size and the last two are. The scale is the average of the normalized counts.

