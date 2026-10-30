# Source and preparation notes

These notes describe the study inputs. To load prepared data with *mutt*, see the
[package README](../../README.md).

## Paper

Study: https://journals.asm.org/doi/full/10.1128/msystems.01182-20

## Files

The published input files T3_SRS_metadata_ms.txt and T4_SRS_final_QF_table_wTaxa.biom were downloaded from: https://github.com/knightlab-analyses/Saliva_quantification_study/tree/master/data

The count table was extracted from the BIOM file using:
```R
library(biomformat)
biom_data <- read_hdf5_biom("T4_SRS_final_QF_table_wTaxa.biom")
Y <- do.call(rbind, biom_data$data)
row.names(Y) <- sapply(biom_data$rows, function(item){item$metadata$Taxon})
saveRDS(Y, "2021_marotz_mSystems_oral_mouthwash.RDS")
```

