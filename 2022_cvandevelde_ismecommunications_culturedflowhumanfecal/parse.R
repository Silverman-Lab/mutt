parse_2022_cvandevelde_ismecommunications_culturedflowhumanfecal <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
            stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                ". Please install them before running this function.")
    }
    if (!is.logical(raw)) {
        stop("raw must be a logical value")
    }
    if (!is.logical(align)) {
        stop("align must be a logical value")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2022_cvandevelde_ismecommunications_culturedflowhumanfecal")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJEB51873_dada2_counts.rds.zip") 
    repro_tax_zip        <- file.path(local, "PRJEB51873_dada2_taxa.rds.zip") 
    scale_16s_zip        <- file.path(local, "VandeVelde2022_scale.csv.zip")
    counts_16s_zip       <- file.path(local, "VandeVelde_2022_16S.csv.zip")
    metadata_16s_zip     <- file.path(local, "VandeVelde_2022_metadata.csv.zip")
    sra_zip              <- file.path(local, "SraRunTable (39).csv.zip")

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    tax_reprocessed <- NA

    # ---- scale and metadata -----
    scale     <- read_zipped_csv(scale_16s_zip) %>% as.data.frame() %>% rename(Sample = SampleID) %>%
                    mutate(log2_FC = ifelse(10^load > 0, log2(10^mean_FC), NA)) %>%
                    rename(log10_FC = load)

    # no way to connect the metadata to the scale and counts because non similar sampleIDs, I suspect it might have to do with plate names but cant figure out.
    # need to email authors about this, but also, there are a ton more samples in the metadata and reprocessed data and in the flow cytometry data that needs to be added.
    # Actually the counts might be in supplementary data, but I dont really know if theyre useable for this -- i dont think so. so need to review study and ask authors about the genomic
    # data connection to our scale. 
    metadata  <- read_zipped_csv(metadata_16s_zip) %>% as.data.frame() %>% rename(Sample = SampleID)
    sra       <- read_zipped_table(sra_zip) %>% rename(Accession = Run) 

    # ------ original counts ------
    counts_original <- read_zipped_csv(counts_16s_zip)

    if (!is.na(counts_original)[1]) {
        original_taxa <- colnames(counts_original)

        # Create taxa mapping data frame
        tax_original <- data.frame(
            Taxa = original_taxa,
            stringsAsFactors = FALSE
        )

        if (!raw) {
            aligned = rename_and_align(counts_original = counts_original, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
            counts_original = aligned$counts_original
        }

        # ------ proportions from counts ------
        proportions_original <- t(counts_original)
        proportions_original <- as.data.frame(t(proportions_original / rowSums(proportions_original, na.rm = TRUE)))
        
    } else {
        proportions_original <- NA
        tax_original <- NA
    }

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
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample_name", align = align, study_name=basename(local))
        counts_reprocessed = aligned$reprocessed
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
        proportions_original = fill_na_zero_numeric(proportions_original)
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    }

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original,
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = proportions_original,
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = tax_original,
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))
}
