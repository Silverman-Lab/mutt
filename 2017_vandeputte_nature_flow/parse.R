parse_2017_vandeputte_nature_flow <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tibble", "tidyverse", "readxl")
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
    library(tibble)
    library(tidyverse)
    library(readxl)

    # ----- Local base directory -----
    local <- file.path("2017_vandeputte_nature_flow")

    # ----- File paths -----
    metadata_zip         <- file.path(local, "Vandeputte_2017_metadata.csv.zip")
    metadata_two_zip     <- file.path(local, "cellcountstotal.csv.zip")
    orig_counts_zip      <- file.path(local, "OTU_nochim.csv.zip")
    orig_tax_rdp_zip     <- file.path(local, "otu_taxonomy_rdp.csv.zip")
    orig_tax_silva_zip   <- file.path(local, "otu_taxonomy_silva.csv.zip")
    orig_prop_zip        <- NA
    repro_counts_rds_zip <- file.path(local, "PRJEB21504_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJEB21504_dada2_taxa.rds.zip")
    sra_zip              <- file.path(local, "SraRunTable (39).csv.zip")

    # ----- Metadata and Scale -----
    # Read scale
    metadata_two <- read_zipped_table(metadata_two_zip, row.names = NULL)  

    # Read metadata
    metadata_df <- read_zipped_table(metadata_zip, row.names = NULL)
    # Remove empty columns and clean up column names
    metadata_df <- metadata_df[, !sapply(metadata_df, function(x) all(is.na(x) | x == ""))]
    # Keep only the first 8 columns which have the actual metadata
    metadata_df <- metadata_df[, 1:8]
    colnames(metadata_df) <- make.unique(colnames(metadata_df))
    
    sra <- read_zipped_table(sra_zip, row.names = NULL) %>% rename(Accession = Run, Sample = Sample_name)

    # Join metadata
    df <- full_join(metadata_df, metadata_two, by = c('Sample' = 'Individual'))
    
    # Create metadata
    metadata <- df %>%
        select(
            Accession = sampleID,
            Sample,
            Cohort = Cohort.x,
            Day,
            `Health status` = Health,
            Enterotype = Enterotype.x
        )

    metadata = merge(metadata, sra, by = "Accession") %>% rename(Sample = Sample.x)

    # Create scale
    scale <- df %>%
        select(
            Sample,
            Day,
            Accession = sampleID,
            `Average cell count (per gram of fresh feces)`,
            `STDEV cell count (per gram of fresh feces)`,
            `Average cell count (per gram of frozen feces)` = `Average cell count (per gram of frozen feces).y`,
            `STDEV cell count (per gram of frozen feces)` = `STDEV cell count (per gram of frozen feces).y`
        ) %>% 
        mutate(log2_FC_cell_g = ifelse(`Average cell count (per gram of fresh feces)`>0, log2(`Average cell count (per gram of fresh feces)`), NA)) %>%
        mutate(log10_FC_cell_g = ifelse(`Average cell count (per gram of frozen feces)`>0, log10(`Average cell count (per gram of frozen feces)`), NA))

    # ----- Original counts from CSV.zip -----
    if (file.exists(orig_counts_zip)) {
    orig_csv <- unzip(orig_counts_zip, list = TRUE)$Name[1]
    orig_con <- unz(orig_counts_zip, orig_csv)
    orig_mat <- read.csv(orig_con, row.names = 1, check.names = FALSE)

    # Read original taxonomy
    orig_tax_rdp_zip <- file.path(local, "otu_taxonomy_rdp.csv.zip")
    orig_tax_rdp_csv <- unzip(orig_tax_rdp_zip, list = TRUE)$Name[1]
    orig_tax_rdp_con <- unz(orig_tax_rdp_zip, orig_tax_rdp_csv)
    orig_tax_rdp <- read.csv(orig_tax_rdp_con, row.names = 1, check.names = FALSE)

    orig_tax_silva_zip <- file.path(local, "otu_taxonomy_silva.csv.zip")
    orig_tax_silva_csv <- unzip(orig_tax_silva_zip, list = TRUE)$Name[1]
    orig_tax_silva_con <- unz(orig_tax_silva_zip, orig_tax_silva_csv)
    orig_tax_silva <- read.csv(orig_tax_silva_con, row.names = 1, check.names = FALSE)

    # Combine taxonomy
    orig_tax <- data.frame(
        Sequence = colnames(orig_mat),
        Kingdom_rdp = orig_tax_rdp$Kingdom,
        Phylum_rdp = orig_tax_rdp$Phylum,
        Class_rdp = orig_tax_rdp$Class,
        Order_rdp = orig_tax_rdp$Order,
        Family_rdp = orig_tax_rdp$Family,
        Genus_rdp = orig_tax_rdp$Genus,
        Kingdom_silva = orig_tax_silva$Kingdom,
        Phylum_silva = orig_tax_silva$Phylum,
        Class_silva = orig_tax_silva$Class,
        Order_silva = orig_tax_silva$Order,
        Family_silva = orig_tax_silva$Family,
        Genus_silva = orig_tax_silva$Genus
    )

    counts_original <- orig_mat
    counts_original$Sequence <- rownames(counts_original)
    counts_original <- counts_original[, c("Sequence", setdiff(names(counts_original), "Sequence"))]
    rownames(counts_original) <- counts_original$Sequence
    } else {
    counts_original <- NA
    }

    if (!raw) {
        align <- rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
        counts_original <- align$counts_original
    }

    if (file.exists(orig_counts_zip)) {
    proportions_original = counts_original
    proportions_original[-1] <- lapply(
      counts_original[-1],
      function(col) col / sum(col)
    )
    } else if (file.exists(orig_prop_zip)) {
    prop_csv <- unzip(orig_prop_zip, list = TRUE)$Name[1]
    prop_con <- unz(orig_prop_zip, prop_csv)
    proportions_original = read.csv(prop_con, row.names = 1, check.names = FALSE) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sequence") %>%
      dplyr::select(Sequence, everything())
    } else {
    proportions_original <- NA
    }

    # --- Original taxonomies ---
    read_taxonomy_zip <- function(zip_path) {
        if (!file.exists(zip_path)) return(NA)
        zip_contents <- unzip(zip_path, list = TRUE)
        csv_name <- zip_contents$Name[1]
        con <- unz(zip_path, csv_name)
        tax_df <- read.csv(con, row.names = 1, check.names = FALSE)
        tax_df$Sequence <- rownames(tax_df)
        tax_df <- tax_df[, c("Sequence", setdiff(names(tax_df), "Sequence"))]
        return(tax_df)
    }

    tax_original_rdp   <- read_taxonomy_zip(orig_tax_rdp_zip)
    tax_original_silva <- read_taxonomy_zip(orig_tax_silva_zip)
    
    # Combine RDP and Silva taxonomies if both are available
    if (!is.na(tax_original_rdp)[1] && !is.na(tax_original_silva)[1]) {
        tax_original <- merge(tax_original_rdp, tax_original_silva, by = "Sequence", suffixes = c("_rdp", "_silva"))
    } else if (!is.na(tax_original_rdp)[1]) {
        tax_original <- tax_original_rdp
    } else if (!is.na(tax_original_silva)[1]) {
        tax_original <- tax_original_silva
    } else {
        tax_original <- NA
    }


    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(repro_counts_rds_zip))) {
        # Create temporary directory for unzipped files
        temp_dir <- tempfile("reprocessed_data")
        dir.create(temp_dir)
        
        # Unzip counts file
        unzip(repro_counts_rds_zip, exdir = temp_dir)
        counts_file <- list.files(temp_dir, pattern = "_counts\\.rds$", full.names = TRUE)[1]
        if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(counts_file))
        
        # Unzip taxonomy file
        unzip(repro_tax_zip, exdir = temp_dir)
        tax_file <- list.files(temp_dir, pattern = "_taxa\\.rds$", full.names = TRUE)[1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))
        
        # Clean up temporary directory
        unlink(temp_dir, recursive = TRUE)
        
        # Convert sequences to lowest rank taxonomy found and update key
        tax_reprocessed = make_taxa_label(tax_reprocessed)
        
        # Convert accessions to sample IDs / Sequences to Taxa
        if (!raw) {
            tryCatch({
                align <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
                counts_reprocessed <- align$counts_reprocessed
                
                # taxa
                matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
                colnames(counts_reprocessed) <- matched_taxa
                
                # Handle duplicate column names by making them unique
                if(any(duplicated(colnames(counts_reprocessed)))) {
                    warning("Duplicate column names found in counts_reprocessed, making them unique")
                    colnames(counts_reprocessed) <- make.unique(colnames(counts_reprocessed))
                }
                
                counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), 
                                                           group = colnames(counts_reprocessed))))
            }, error = function(e) {
                warning("Error processing reprocessed data: ", e$message)
            })
        }
        
        # Calculate proportions
        if (!is.na(counts_reprocessed)[1]) {
            tryCatch({
                proportions_reprocessed = counts_reprocessed
                proportions_reprocessed[-1] <- lapply(
                    counts_reprocessed[-1],
                    function(col) col / sum(col)
                )
            }, error = function(e) {
                warning("Error calculating reprocessed proportions: ", e$message)
            })
        }
    } else {
        counts_reprocessed <- NA
        tax_reprocessed <- NA
        proportions_reprocessed <- NA
    }

    if (!raw) {
        tryCatch({
            if (!is.na(counts_original)[1]) {
                counts_original = fill_na_zero_numeric(counts_original)
            }
            if (!is.na(proportions_original)[1]) {
                proportions_original = fill_na_zero_numeric(proportions_original)
            }
            if (!is.na(counts_reprocessed)[1]) {
                counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
            }
            if (!is.na(proportions_reprocessed)[1]) {
                proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
            }
        }, error = function(e) {
            warning("Error filling NA values: ", e$message)
        })
    }

    # ----- Return all -----
    return(list(
      counts      = list(
        original    = counts_original,
        reprocessed = counts_reprocessed
      ),
      proportions = list(
        original    = proportions_original,
        reprocessed = proportions_reprocessed
      ),
      tax         = list(
        original = tax_original,
        reprocessed = tax_reprocessed
      ),
      scale       = scale,
      metadata    = metadata
    ))
}
