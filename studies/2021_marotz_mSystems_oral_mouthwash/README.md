## Paper

This is from the paper: https://journals.asm.org/doi/full/10.1128/msystems.01182-20

## Files

The raw files T3_SRS_metadata_ms.txt and T4_SRS_final_QF_table_wTaxa.biom were downloaded from: https://github.com/knightlab-analyses/Saliva_quantification_study/tree/master/data

The biom file was processed as
```R
library(biomformat)
biom_data <- read_hdf5_biom("T4_SRS_final_QF_table_wTaxa.biom")
Y <- do.call(rbind, biom_data$data)
row.names(Y) <- sapply(biom_data$rows, function(item){item$metadata$Taxon})
saveRDS(Y, "2021_marotz_mSystems_oral_mouthwash.RDS")
```

