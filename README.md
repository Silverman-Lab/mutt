# data_repository
data respository for the Silverman lab

Before download, **make sure gitlfs is installed**. 

Use Git LFS for large files: https://docs.github.com/en/repositories/working-with-files/managing-large-files/configuring-git-large-file-storage

According to Maxwell, metagenomic data is huge. So lets keep this to count tables for the moment (already processed data). 

# Notes
- A growing list of studies and relevant metadata is being maintained here: https://docs.google.com/spreadsheets/d/13b4Toscse0MjyAGYt1zfWoPxSRpyuvENHVGKBLYwAAw/edit?usp=sharing
- An extra column in a count table that contains non-meaningful information (e.g., all values are -1) and is not present in the taxa table indicates that the corresponding taxa are unassigned.

# Format to Maintain
- Each dataset gets its own directory with an informative name (all lower-case), "year_name_journal_keyword". Try to be concise, but not too vague (don't use date ranges e.g., 2020-2024 but just choose a single date for the study). 

- Datasets should be compressed (and again, stored using Git LFS: https://docs.github.com/en/repositories/working-with-files/managing-large-files/configuring-git-large-file-storage). This can be done by running ./zip-push-gitlfs.sh from the terminal within the data repository main folder.

## Parsed Data Structure for `parse.R` scripts. 
parse.R should have a single function named `parse_[name of directory]` (all lowercase) which returns a list object with the following elements. That function should not require arguments but they can be optional. Ideally, parse scripts use nothing other than base R or tidyverse functions to minimize dependencies and errors if certain libraries are not installed. 

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
`parse.R` scripts should specify file paths relative to the root `data_repository` directory. 

### Helper Functions:
`repoparser.py` a parser check function to check the structure and execution of each individual parser. (Needs work)
`obtainpublicationinfo_pmid.py` using a list of PMIDs, this functionality can be integrated into each parser to obtain the manuscript information from NCBI (Script works, but python and each parse script is in R.)
`zip-push-gitlfs.sh` run from terminal in the repository directory when you are ready to push and it will compress your files with .zip and upload with gitlfs

### Wrapper functions for MicrobialScaleRepository package:
`microbialscalerepo.R` function to call parse scripts (with selection of individual studies) and optionally store in .Rdata object
`microbialscalerepo.py` function to call parse scripts (with selection of individual studies) and optionally store in .pkl object

## For now, functionality:

```r

# Source the files from the directory -- No R package just yet
source("data_repository/microbialscalerepo.R") 
source("data_repository/helpers.R")

# Choose whichever dataset you want and supply like this or as named vector, or just supply a vector of the repo directories:
study_parsers <- c(
   Vandeputte2021 = "2021_vandeputte_naturecommunications_flow_timeseries",
   CvandeVelde2022 = "2022_cvandevelde_ismecommunications_culturedflowhumanfecal",
   Vandeputte2017 = "2017_vandeputte_nature_flow",
   Pereira2023 = "2023_pereira_nature_nervous",
   Krawczyk2022 = "2022_krawczyk_microbiome_tickgeographicaldistributionqpcr",
   Liao2021 = "2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct",
   Stammler2016 = "2016_stammler_microbiome_micehuman",
   Dreier2022 = "2022_dreier_bmcmicrobiology_cheeseqpcr",
   GALAXY = "2024_nishijima_cell_galaxy",
   MetaCardis = "2022_fromentin_naturemedicine_metacardissubset",
   Marotz2021 = "2021_marotz_mSystems_oral_mouthwash",
   Morton2019 = "2019_morton_naturecommunications_songbird_oral",
   Vieira_Silva2019 = "2019_vieirasilva_naturemicrobiology_pscibd",
   Contijoch2019 = "2019_contijoch_elife_multispeciesqPCRshotgunandamplicon",
   Alessandri2024 = "2024_alessandri_microbbiotechnology_pcosvaginalmicrobiota",
   Maghini2023 = "2023_maghini_naturebiotechnology_samplemesurement",
   Garcia_Martinez2024 = "2024_garciamartinez_bmcmicrobiology_ckdanddysbiosiswithserum",
   Sternes2024 = "2024_sternes_frontmicrobiol_IBDppiqPCR",
   Rao2021 = "2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification",
   Tettamanti_Boshier2020 = "2020_tettamantiboshier_msystems_vaginaltimeseries",
   Kruger2024 = "2024_kruger_scientificreports_ddpcrhealthysubjects",
   Liu2017 = "2017_liu_mbio_penilehivqPCR",
   Fu2023 = "2023_fu_imeta_wasterwater_pathogens",
   Jin2022 = "2022_jin_natureComm_technicalReplicates",
   Zaramela2022 = "2022_zaramela_msystems_synDNA",
   Feng2023 = "2023_feng_imetawiley_chickensegment",
   Reese2022 = "2021_reese_cell_chimpanzee",
   Barlow2020 = "2020_barlow_naturecommunications_miceGI",
   Morton2019 = "2019_morton_naturecommunications_songbird_oral",
   Prochazkova2024 = "2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal",
   Zemb2020 = "2020_zemb_microOpen_spike",
   Jin2024 = "2024_jin_pnas_semen",
   Galazzo2020 = "2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy",
   Lin2019 = "2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein",
   Suriano2022 = "2022_suriano_aps_micefecal",
   Thiruppathy2025 = "2025_thiruppathy_microbiome_relicDNAflow",
   Wagner2025 = "2025_wagner_frontiersinmiccrobiology_flowpiglets",
   Kallastu2023 = "2023_kallastu_research_foodscience_food"
)

# Run repo function
studies <- microbialscalerepo(
  studies = study_parsers, # If not supplied, defaults to all
  base_directory = "data_repository/", # This is default, but you should change to wherever your local download is, for now.
  rawdata = FALSE, # Dont change this because its the un-reformatted original data non cleaned. If TRUE, returns unformatted original data
  align_samples = FALSE, # If TRUE, this will align your matrices to the scale dataframe so all sample data is aligned (If it can be)
  save_to = "datasetsfromrepo.RData", # OPTIONAL, save RData object of all the studies you chose.
  verbose = TRUE # Display structure of datasets returned if TRUE
)
```

