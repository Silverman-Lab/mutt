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
parse.py should have a single function named `parse_[name of directory]` which returns a list object with the following elements. That function should not require arguments but they can be optional. Ideally, parse scripts use nothing other than base R or tidyverse functions to minimize dependencies and errors if certain libraries are not installed. 

- `counts` integer valued count matrix (not data.frame) that is (D x N) and has row and column names which are sampleIDs and sequenceIDs (e.g., taxaIDs) respectively. 
- `proportions` real-valued valued count matrix (not data.frame) that is (D x N) and has row and column names which are sampleIDs and sequenceIDs (e.g., taxaIDs) respectively. only include if counts are not available
- `scale` positive-valued matrix likely of dimension 1 x N but other formats may need to be allowed. Column names should be sampleIDs
- `metadata` (optional but often required) N x Q data.frame rownames should be sampleIDs. 
- `tax` (optional) D x ?, character-valued data.frame with sequenceIDs as rownames and each column labeled in a meaningful way. For microbiome data these labels should be limited to c("kingdom", "phylum", "class", "order", "genus", "species", "strain")
- `sequences` (optional) D x ?, character valued data.frame or matrix linking sequenceIDs (rownames) to the actual raw sequence (e.g., the 16S sequence of a particular taxon)
- `phylo` (optional) phylogenetic tree stored in reasonable format (let me know if any repos have phylogenetic trees in them and I will figure out a good standard format)

### Note on file paths
`parse.R` scripts should specify file paths relative to the root `data_repository` directory. See my example in `2024_nishijima_cell_galaxy/parse.R`. 
