datasets <- list(
  "2021_vandeputte_naturecommunications_flow_timeseries" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = NULL,
    ID             = "Sample_ID",
    PMID           = "34795283"
  ),
  "2022_cvandevelde_ismecommunications_culturedflowhumanfecal" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = "Prefix",
    ID             = "Sample",
    PMID           = "37938658"
  ),
  "2017_vandeputte_nature_flow" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = "Status",
    ID             = "Sample",
    PMID           = "29143816"
  ),
  "2023_pereira_nature_nervous" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("condition", "time_hours", "medium", "replicate"),
    ID             = "Sample",
    PMID           = "39572788"
  ),
  "2022_krawczyk_microbiome_tickgeographicaldistributionqpcr" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Tick",
    sampletype     = "Whole Organism",
    covariates     = c("Life stage", "Location"),
    ID             = "Sample_name",
    PMID           = "35927748"
  ),
  "2021_liao_scientificdata_longitudinalmicrobiomeqpcr_allohct" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("VanA", "Consistency", "Pool", "Timepoint"),
    ID             = "SampleID",
    PMID           = "33654104"
  ),
  "2016_stammler_microbiome_micehuman" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Mice",
    sampletype     = "Fecal",
    covariates     = c("Treatment", "Dilution", "Background"),
    ID             = "SampleID",
    PMID           = "27329048"
  ),
  "2022_dreier_bmcmicrobiology_cheeseqpcr" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Cheese",
    sampletype     = "Cheese",
    covariates     = NULL,
    ID             = "Sample_name",
    PMID           = "35130830"
  ),
  "2024_nishijima_cell_galaxy" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = "cohort",
    ID             = "ID",
    PMID           = "39541968"
  ),
  "2022_fromentin_naturemedicine_metacardissubset" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("PatientGroup", "SMOKE", "Gender", "CENTER", "AGE", "cohort", "pa_work_2cl"),
    ID             = "Sample",
    PMID           = c("35177860", "39541968", "34880489")
  ),
  "2021_marotz_mSystems_oral_mouthwash" = list(
    sequencingtype = "16S rRNA",
    loadtype       = c("Flow Cytometry", "qPCR"),
    organismtype   = "Human",
    sampletype     = "Saliva",
    covariates     = c("Treatment_code", "sex", "host_age", "mouthwash_regularly", "processing"),
    ID             = "SampleID",
    PMID           = "33594005"
  ),
  "2019_vieirasilva_naturemicrobiology_pscibd" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("Diagnosis", "Age", "Gender", "BMI"),
    ID             = "ID",
    PMID           = "31209308"
  ),
  "2019_contijoch_elife_multispeciesqPCRshotgunandamplicon" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype       = "qPCR",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("experimental_timepoint", "sample_mass"),
    ID             = "Sample Name",
    PMID           = "30666957"
  ),
  "2024_tunsakul_peerj_aerobicvsanaerobicinhealthyvsobesity" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("environment"),
    ID             = "Sample",
    PMID           = "38650647"
  ),
  "2024_alessandri_microbbiotechnology_pcosvaginalmicrobiota" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = c("Fecal", "Vaginal"),
    covariates     = c("treatmentfromlibraryname", "env_broad_scale"),
    ID             = "Sample",
    PMID           = "39364592"
  ),
  "2023_maghini_naturebiotechnology_samplemesurement" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype       = "qPCR",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("Condition", "Donor"),
    ID             = "ID",
    PMID           = "37106038"
  ),
  "2024_garciamartinez_bmcmicrobiology_ckdanddysbiosiswithserum" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("population"),
    ID             = "Sample",
    PMID           = "39462312"
  ),
  "2024_sternes_frontmicrobiol_IBDppiqPCR" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Human",
    sampletype     = "Gastrointestinal mucosa biopsy",
    covariates     = c("AGE", "Sex", "BMI", "PPI", "NC_total_score"),
    ID             = "Sample_name",
    PMID           = "39469457"
  ),
  "2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification" = list(
    sequencingtype = "16S rRNA",
    loadtype       = c("Flow Cytometry", "qPCR", "MK_spike"),
    organismtype   = c("Human", "Mouse", "Mock"),
    sampletype     = "Fecal",
    covariates     = c("sample_type"),
    ID             = "Sample_name",
    PMID           = "33627867"
  ),
  "2020_tettamantiboshier_msystems_vaginaltimeseries" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Human",
    sampletype     = "Vaginal muscosa",
    covariates     = c("Hours_In_Study"),
    ID             = "Sample_ID",
    PMID           = "32265316"
  ),
  "2024_kruger_scientificreports_ddpcrhealthysubjects" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "ddPCR",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("BristolStoolScale_highest", "Timepoint", "pH"),
    ID             = "SampleID",
    PMID           = "39427011"
  ),
  "2017_liu_mbio_penilehivqPCR" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Human",
    sampletype     = "Coronal sculcus",
    covariates     = c("HIV Status"),
    ID             = "Sample_name",
    PMID           = "28743816"
  ),
  "2022_jin_natureComm_technicalReplicates" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "ddPCR",
    organismtype   = "Mice",
    sampletype     = "Cecum content",
    covariates     = c("diet", "location", "sample_type", "subject"),
    ID             = "Sample",
    PMID           = "35194029"
  ),
  "2022_zaramela_msystems_synDNA" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Saliva",
    covariates     = c("gender", "Pool"),
    ID             = "ID",
    PMID           = "36317886"
  ),
  "2023_feng_imetawiley_chickensegment" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Chicken",
    sampletype     = "Gut segment",
    covariates     = c("Segment", "Date"),
    ID             = "Sample",
    PMID           = "38868437"
  ),
  "2021_reese_cell_chimpanzee" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Chimpanzee",
    sampletype     = "Fecal",
    covariates     = c("Age group", "Community", "Sex", "YearQuarter"),
    ID             = "Sample_name",
    PMID           = "33232664"
  ),
  "2020_barlow_naturecommunications_miceGI" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "ddPCR",
    organismtype   = "Mice",
    sampletype     = "Fecal",
    covariates     = c("Diet", "Site", "Cage.y", "mouse.y"),
    ID             = "Sample",
    PMID           = "32444602"
  ),
  "2019_morton_naturecommunications_songbird_oral" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype       = c("qPCR", "Flow Cytometry"),
    organismtype   = "Human",
    sampletype     = "Saliva",
    covariates     = c("processing", "sex", "Host_age", "mouthwash_regularly", "treatment"),
    ID             = "Accession",
    PMID           = "31222023"
  ),
  "2024_prochazkova_naturemicrobiology_longitudinalhealthyflowfecal" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("Bristol_scale_Mean", "Faecal_pH_Mean", "Stool_freq_Mean", "Stool_moisture_Mean"),
    ID             = "ID",
    PMID           = "39604623"
  ),
  "2020_zemb_microOpen_spike" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Pig",
    sampletype     = "Fecal",
    covariates     = c("mgfeces", "added_Coli"),
    ID             = "Sample.name",
    PMID           = "31927795"
  ),
  "2024_jin_pnas_semen" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "qPCR",
    organismtype   = "Human",
    sampletype     = "Semen",
    covariates     = c("Fertility", "Diagnosis", "Group", "Sperm_concentration", "Sperm_motility", "Age"),
    ID             = "Sample",
    PMID           = "39359400"
  ),
  "2020_galazzo_frontiersincellularandinfectionmicrobiology_flowqPCRddPCRhealthy" = list(
    sequencingtype = "16S rRNA",
    loadtype       = c("qPCR", "ddPCR", "Flow Cytometry"),
    organismtype   = "Human",
    sampletype     = "Fecal",
    covariates     = c("sample_treatment", "host_subject_id"),
    ID             = "Sample_name",
    PMID           = "32850498"
  ),
  "2019_lin_applenvironmicrobiol_16s18smarineecologyflowandspikein" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Marine",
    sampletype     = "Seawater",
    covariates     = c("Station", "Line", "Filtered Seawater Vol [L]", "SurfChl [mg m-3]"),
    ID             = "SampleID",
    PMID           = "30552195"
  ),
  "2022_suriano_aps_micefecal" = list(
    sequencingtype = "16S rRNA",
    loadtype       = "Flow Cytometry",
    organismtype   = "Mice",
    sampletype     = "Fecal",
    covariates     = c("Day"),
    ID             = "Sample_name",
    PMID           = "36516223"
  ),
  "2025_thiruppathy_microbiome_relicDNAflow" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype       = "Flow Cytometry",
    organismtype   = "Human",
    sampletype     = "Skin",
    covariates     = c("PME_treated", "BodySite_ID", "Skin_Type", "Face_wash", "Exercise_per_week", "Face_sunscreen", "Body_shower_fq_by_day", "Age", "Sex", "Ethnicity", "Swab_Dimensions_sqcm"),
    ID             = "Sample",
    PMID           = "40038838"
  ),
  "2023_kallastu_research_foodscience_food" = list(
    sequencingtype = "16S rRNA",
    loadtype       = c("Flow Cytometry", "qPCR", "CFU"),
    organismtype   = "Food",
    sampletype     = "Food",
    covariates     = c("pma_treatment", "Organism"),
    ID             = "Sample",
    PMID           = "36691592"
  )
)
