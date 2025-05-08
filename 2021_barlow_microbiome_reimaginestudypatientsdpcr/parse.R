parse_2021_barlow_microbiomee_reimaginestudypatientsdpcr <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop(
        "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
        ". Please install them before running this function."
        )
    }
    if (!is.logical(raw) || length(raw) != 1) {
        stop("`raw` must be a single logical value (TRUE or FALSE)")
    }
    if (!is.logical(align) || length(align) != 1) {
        stop("`align` must be a single logical value (TRUE or FALSE)")
    }

    library(tidyverse)
    library(readxl)
    library(stringr)
    library(readr)
    
    # ------- local path -----------------
    localPath <- file.path("2021_barlow_microbiomee_reimaginestudypatientsdpcr")

    # ------- file path -----------------
    repro_counts_rds_zip<- file.path(localPath, "PRJNA674353_dada2_counts.rds.zip")
    repro_tax_zip       <- file.path(localPath, "PRJNA674353_dada2_taxa.rds.zip")
    metadata_sra_zip    <- file.path(local, "SraRunTable (33).csv.zip")
    scale_zip           <- file.path(local, "")

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA
    scale = NA
    metadata = NA
    
    # ------- scale ---------------------
    # ------- metadata ------------------
    sra <- read_zipped_table(metadata_sra_zip) %>% rename(Accession = Run)
    
    # ------- original data -------------
    ########## NEEDS SOMEONE TO WORK ON#########
    # Have to run the jupyter notebooks to get the original data and the scale etc.  
    # ------- counts --------------------
    # ------- proportions ---------------
    # ------- taxa ----------------------

    # ----- Reprocessed counts from RDS ZIP -----
    temp_rds <- tempfile(fileext = ".rds")
    unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

    rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
    if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
    counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

    # ----- Taxonomy reprocessed -----
    temp_tax <- tempfile(fileext = ".rds")
    unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

    tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
    if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

    
    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
        align <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "SampleID", align = align, study_name = basename(local))
        counts_reprocessed <- align$reprocessed
    }

    # taxa
    if (!raw) {
        matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
        colnames(counts_reprocessed) <- matched_taxa
        counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
    }

    # proportions reprocessed
    proportions_reprocessed = counts_reprocessed
    proportions_reprocessed[-1] <- lapply(
        counts_reprocessed[-1],
        function(col) col / sum(col)
    )

    if (!raw) {
      counts_original = fill_na_zero_numeric(counts_original)
      counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
      proportions_original = fill_na_zero_numeric(proportions_original)
      proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }
    
    return(list(
        scale=scale, 
        metadata=metadata, 
        counts=list(
                    original = counts_original, 
                    reprocessed = counts_reprocessed
                ),
        proportions=list(
                    original = proportions_original,
                    reprocessed = proportions_reprocessed
                ),
        tax=list(
                    original = taxa_original,
                    reprocessed = tax_reprocessed
                )
        )
    )

}