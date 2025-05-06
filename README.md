# data_repository
data respository for the Silverman lab

Use Git LFS for large files: https://docs.github.com/en/repositories/working-with-files/managing-large-files/configuring-git-large-file-storage

According to Maxwell, metagenomic data is huge. So lets keep this to count tables for the moment (already processed data). 

# Notes
- A growing list of studies and relevant metadata is being maintained here: https://docs.google.com/spreadsheets/d/13b4Toscse0MjyAGYt1zfWoPxSRpyuvENHVGKBLYwAAw/edit?usp=sharing
- An extra column in a count table that contains non-meaningful information (e.g., all values are -1) and is not present in the taxa table indicates that the corresponding taxa are unassigned.

# Format to Maintain
- Each dataset gets its own directory with an informative name (all lower-case), "year_name_journal_keyword". Try to be concise (don't use date ranges e.g., 2020-2024 but just choose a single date for the study). 

- datasets should be compressed (and again, stored using Git LFS: https://docs.github.com/en/repositories/working-with-files/managing-large-files/configuring-git-large-file-storage)

## Parsed Data Structure for `parse.R` scripts. 
parse.R should have a single function named `parse_[name of directory]` which returns a list object with the following elements. That function should not require arguments but they can be optional. Ideally, parse scripts use nothing other than base R or tidyverse functions to minimize dependencies and errors if certain libraries are not installed. 

- `counts` integer valued count matrix (not data.frame) that is (N x D). sampleIDs (rows) and sequenceIDs (columns) (e.g., taxaIDs) respectively. Should contain a column key with sampleIDs linking to proportions, scale, and metadata.
- `proportions` real-valued valued count matrix (not data.frame) that is (N x D) and has row and column names which are sampleIDs and sequenceIDs (e.g., taxaIDs) respectively. Should contain a column key with sampleIDs linking to counts, scale, and metadata. 
- `scale` positive-valued matrix likely of dimension N x 1 but other formats may need to be allowed due to mean and sd or multiple techniques measuring total scale. Should contain a column key with sampleIDs linking to counts, proportions, and metadata. 
- `metadata` (optional but often required) N x Q data.frame. Should contain a column key with sampleIDs linking to counts and proportions. 
- `tax` (optional) D x ?, character-valued data.frame with sequenceIDs as rownames and each column labeled in a meaningful way. For microbiome data these labels should be limited to c("Kingdom", "Phylum", "Class", "Order", "Genus", "Species", "Strain"). "Taxa" is the lowest identified taxonomy classified specified by prefix and then the classified taxa, if unclassified by lowest taxonomy resolution then prefixed with uc_ and then taxonomic level prefix ex. for phylum: uc_p_[taxa classification name] "Sequence" column links sequenceIDs (ASV, OTUs, classified taxa) to the actual raw sequence (e.g., the 16S sequence of a particular taxon).
- `phylo` (optional) phylogenetic tree stored in reasonable format (let me know if any repos have phylogenetic trees in them and I will figure out a good standard format)

```r
# ----- Example Shotgun Metagenomics Study -----
return(list(
    counts = list(
      original = counts_original,
      reprocessed = list(
          mOTU3 = mOTU3_counts,
          MetaPhlAn4 = MetaPhlAn4_counts
      )
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = list(
          mOTU3 = mOTU3_proportions,
          MetaPhlAn4 = MetaPhlAn4_proportions
      )
    ),
    tax = list(
      original = tax_original,
      reprocessed = list(
          mOTU3 = mOTU3_tax,
          MetaPhlAn4 = MetaPhlAn4_tax
      )
    ),
    scale = scale,
    metadata = metadata,
    phylo = NA
))

# ---- Example amplicon study ----
return(list(
    counts = list(
        original = counts_original,
        reprocessed = counts_reprocessed
    ),
    tax = list(
        original = tax_original,
        reprocessed = tax_reprocessed
    ),
    proportions = list(
        original = proportions_original,
        reprocessed = proportions_reprocessed
    ),
    metadata = metadata,
    scale = scale,
    phylo = NA
))
```
### Note on file paths
`parse.R` scripts should specify file paths relative to the root `data_repository` directory. See my example in `2024_nishijima_cell_galaxy/parse.R`. 
