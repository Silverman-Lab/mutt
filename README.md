# Microbial Scale Data Repository for the Silverman lab
<br>
Authors: Maxwell Konnaris, Justin Silverman MD PhD <br>
Affiliation: Pennsylvania State University
<br><br>

Before download, **make sure gitlfs is installed**. 

## Notes
- A growing list of studies and relevant metadata is being maintained here: https://docs.google.com/spreadsheets/d/13b4Toscse0MjyAGYt1zfWoPxSRpyuvENHVGKBLYwAAw/edit?usp=sharing
- Use Git LFS for large files: https://docs.github.com/en/repositories/working-with-files/managing-large-files/configuring-git-large-file-storage
- Shotgun metagenomic data is huge, therefore, this is limited to processed data tables for the moment. 

## Format to Maintain
- Each dataset gets its own directory with an informative name (all lower-case), "year_name_journal_keyword". Try to be concise, but not too vague (don't use date ranges e.g., 2020-2024 but just choose a single date for the study). 

- Datasets should be compressed before uploading (and again, stored using Git LFS: https://docs.github.com/en/repositories/working-with-files/managing-large-files/configuring-git-large-file-storage). This can be done by running ./zip-push-gitlfs.sh from the terminal within the data repository main folder while pushing.

## Parsed Data Structure for `parse.R` scripts. 
parse.R should have a single function named `parse_[name of directory]` (all lowercase) which returns a list object with the following elements. That function should not require arguments but they can be optional. Ideally, parse scripts use nothing other than base R or tidyverse functions to minimize dependencies and errors if certain libraries are not installed. 

- `counts` integer valued count matrix (not data.frame) that is (N x D). sampleIDs (rows) and sequenceIDs (columns) (e.g., taxaIDs) respectively. Should contain a column key with sampleIDs linking to proportions, scale, and metadata.
- `proportions` real-valued valued count matrix (not data.frame) that is (N x D) and has row and column names which are sampleIDs and sequenceIDs (e.g., taxaIDs) respectively. Should contain a column key with sampleIDs linking to counts, scale, and metadata. 
- `scale` positive-valued matrix likely of dimension N x 1 but other formats may need to be allowed due to mean and sd or multiple techniques measuring total scale. Should contain a column key with sampleIDs linking to counts, proportions, and metadata. 
- `metadata` (optional but often required) N x Q data.frame. Should contain a column key with sampleIDs linking to counts and proportions. 
- `tax` (optional) D x ?, character-valued data.frame with sequenceIDs as rownames and each column labeled in a meaningful way. For microbiome data these labels should be limited to c("Kingdom", "Phylum", "Class", "Order", "Genus", "Species", "Strain"). "Taxa" is the lowest identified taxonomy classified specified by prefix and then the classified taxa, if unclassified by lowest taxonomy resolution then prefixed with uc_ and then taxonomic level prefix ex. for phylum: uc_p_[taxa classification name] "Sequence" column should additionally be the tax rownames and links sequenceIDs (ASV, OTUs, classified taxa) to the actual raw sequence (e.g., the 16S sequence of a particular taxon).
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
## Note on file paths
`parse.R` scripts should specify file paths relative to the root `data_repository` directory. 

## Helper Scripts:
- `obtainpublicationinfo_pmid.py` using a list of PMIDs, this functionality can be integrated into each parser to obtain the manuscript information from NCBI (Script works, but python and each parse script is in R.)
- `zip-push-gitlfs.sh` run from terminal in the repository directory when you are ready to push and it will compress your files with .zip and upload with gitlfs
- `downloadrepoaspkl.py` convert the repo to a .pkl object
- `loadRDataintopython.ipynb` example using rpy2 to load the RData object returned from `microbialscalerepo.R` into python

## Wrapper functions for MicrobialScaleRepository package:
- `microbialscalerepo.R` function to call parse scripts (with selection of individual studies) and optionally store in .Rdata object

# Functionality (for now):

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
    Vieira_Silva2019 = "2019_vieirasilva_naturemicrobiology_pscibd",
    Contijoch2019 = "2019_contijoch_elife_multispeciesqPCRshotgunandamplicon",
    Tunsakul2024 = "2024_tunsakul_peerj_aerobicvsanaerobicinhealthyvsobesity",
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
    #Barlow2021 = "2021_barlow_microbiome_reimaginestudypatientsddpcr",
    Morton2019 = "2019_morton_naturecommunications_songbird_oral",
    Prochazkova2024 = "2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal",
    Zemb2020 = "2020_zemb_microOpen_spike",
    Jin2024 = "2024_jin_pnas_semen",
    Galazzo2020 = "2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy",
    Lin2019 = "2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein",
    Suriano2022 = "2022_suriano_aps_micefecal",
    Thiruppathy2025 = "2025_thiruppathy_microbiome_relicDNAflow",
    Wagner2025 = "2025_wagner_frontiersinmiccrobiology_flowpiglets",
    #Rolling2021 = "2021_rolling_naturemicrobiology_allohct_ITS_16s",
    Kallastu2023 = "2023_kallastu_research_foodscience_food"
)

# Run repo function
repo <- microbialscalerepo(
  studies = study_parsers, # If not supplied, defaults to all
  base_directory = "data_repository/", # This is default, but you should change to wherever your local download is, for now.
  rawdata = FALSE, # Dont change this because its the un-reformatted original data non cleaned. If TRUE, returns unformatted original data
  align_samples = FALSE, # If TRUE, this will align your matrices to the scale dataframe so all sample data is aligned (If it can be)
  save_to = "datasetsfromrepo.RData", # OPTIONAL, save RData object of all the studies you chose.
  verbose = TRUE # Display structure of datasets returned if TRUE
)
```

This will supply you with the three saved objects:
- *.RData object containing all of the parsed data formatted from each study
- *_validation.RData containing the validated structure of the parsed data from each study contained in the .RData object
- *_validation_summary.txt containing the skeleton summarized structure of what is within the .RData object and if verbose = TRUE, this is printed to console. See example below:

```r
Study: Vandeputte2017
Structure:
  counts: LIST={original = LIST={rdp19 = DATAFRAME(135 x 193), rdp16 = DATAFRAME(135 x 140)}, nishijima2024 = DATAFRAME(95 x 166), reprocessed = LIST={rdp19 = DATAFRAME(135 x 334), rdp16 = DATAFRAME(135 x 222)}}
  proportions: LIST={original = LIST={rdp19 = DATAFRAME(135 x 193), rdp16 = DATAFRAME(135 x 140)}, nishijima2024 = DATAFRAME(95 x 166), reprocessed = LIST={rdp19 = DATAFRAME(135 x 334), rdp16 = DATAFRAME(135 x 222)}}
  tax: LIST={original = LIST={rdp19 = DATAFRAME(2041 x 14), rdp16 = DATAFRAME(2041 x 7)}, nishijima2024 = DATAFRAME(166 x 1), reprocessed = LIST={rdp19 = DATAFRAME(5173 x 7), rdp16 = DATAFRAME(5173 x 7)}}
  scale: DATAFRAME(321 x 11)
  metadata: DATAFRAME(321 x 56)
  studydemographics: LIST={sequencingtype = VECTOR(length=1), loadtype = VECTOR(length=1), organismtype = VECTOR(length=1), sampletype = VECTOR(length=1), covariates = VECTOR(length=1), ID = VECTOR(length=1), PMID = VECTOR(length=1)}
```

Additionally, there are many callable helper functions to scrape additional data about each study in the returned .RData object:

```r
# WORK IN PROGRESS - Works for obtainpublicationinfo_pmid.py

# library(reticulate) is required for python packages
annotate_studies_with_pubmed <- function(repo,
                                         script_path, # Supply any script (R or Python) that works off the structure of the repo such as obtainpublicationinfo_pmid.py
                                         email         = NA_character_,
                                         api_key       = NA_character_,
                                         overwrite     = TRUE,
                                         csv_copy_to   = "publication_data.csv") 
```

`annotate_studies_with_pubmed()` will append data to the repo for each study as a list
The data appended is determined by the returned object from the script. The intention is modularity for users outside of provided scripts in the repo readily available. The example, `obtainpublicationinfo_pmid.py` will scrape information about the published journal article from NCBI and append to $studyinfo. See below:

```r

```

