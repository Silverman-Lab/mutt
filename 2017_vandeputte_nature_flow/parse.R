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
    metadata_two <- read_zipped_table(metadata_two_zip, row.names = NULL) %>% rename(Sample = Individual)

    # Read metadata
    metadata_df <- read_zipped_table(metadata_zip, row.names = NULL)
    metadata_df <- metadata_df[, !sapply(metadata_df, function(x) all(is.na(x) | x == ""))]
    metadata_df <- metadata_df[, 1:8]
    colnames(metadata_df) <- make.unique(colnames(metadata_df))
    
    sra <- read_zipped_table(sra_zip, row.names = NULL) %>% rename(Accession = Run, Sample = Sample_name)

    # Join metadata
    df <- full_join(metadata_df, metadata_two, by = 'Sample')
    
    # Create metadata
    metadata <- df %>%
        select(
            Accession = sampleID,
            Sample,
            Cohort = Cohort.x,
            Day,
            `Health status` = Health,
            Enterotype = Enterotype.x
        ) %>%
        as.data.frame()  

    metadata <- merge(metadata, sra, by = "Accession", suffixes = c(".x", ".y")) %>% rename(Sample = Sample.x)

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
        mutate(log10_FC_cell_g = ifelse(`Average cell count (per gram of frozen feces)`>0, log10(`Average cell count (per gram of frozen feces)`), NA)) %>%
        as.data.frame()  

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
            Kingdom= orig_tax_rdp$Kingdom,
            Phylum = orig_tax_rdp$Phylum,
            Class = orig_tax_rdp$Class,
            Order = orig_tax_rdp$Order,
            Family = orig_tax_rdp$Family,
            Genus = orig_tax_rdp$Genus,
            Kingdom_silva = orig_tax_silva$Kingdom,
            Phylum_silva = orig_tax_silva$Phylum,
            Class_silva = orig_tax_silva$Class,
            Order_silva = orig_tax_silva$Order,
            Family_silva = orig_tax_silva$Family,
            Genus_silva = orig_tax_silva$Genus
        )

        tax_original <- make_taxa_label(orig_tax)
        tax_original <- tax_original %>% rename(Kingdom_rdp = Kingdom, Phylum_rdp = Phylum, Class_rdp = Class, Order_rdp = Order, Family_rdp = Family, Genus_rdp = Genus)
        rownames(tax_original) <- tax_original$Sequence
        counts_original <- orig_mat

        if (!raw) {
            aligned <- rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
            counts_original <- aligned$counts_original
            matched_taxa <- tax_original$Taxa[match(colnames(counts_original), rownames(tax_original))]
            colnames(counts_original) <- matched_taxa
            counts_original <- collapse_duplicate_columns_exact(counts_original)
            original_names <- colnames(counts_original)
            counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
        }

        # Ensure counts are numeric after alignment
        counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original))
        proportions_original <- sweep(counts_original, MARGIN = 1, STATS = rowSums(counts_original), FUN = "/")

    } else if (file.exists(orig_prop_zip)) {
        counts_original <- NA
        prop_csv <- unzip(orig_prop_zip, list = TRUE)$Name[1]
        prop_con <- unz(orig_prop_zip, prop_csv)
        proportions_original = read.csv(prop_con, row.names = 1, check.names = FALSE) %>%
            as.data.frame() %>%
            tibble::rownames_to_column("Sequence") %>%
            dplyr::select(Sequence, everything()) %>% t() %>% as.data.frame()

            proportions_original <- as.data.frame(lapply(proportions_original, as.numeric), row.names = rownames(proportions_original))

    } else {
        counts_original <- NA
        proportions_original <- NA
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
        
        # Ensure counts are numeric
        counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed))
        
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
                aligned <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
                counts_reprocessed <- aligned$reprocessed
                
                # taxa
                matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
                colnames(counts_reprocessed) <- matched_taxa
                counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
                original_names <- colnames(counts_reprocessed)
                counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
                
                # Calculate proportions using sweep
                proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
            }, error = function(e) {
                warning("Error in processing reprocessed data: ", e$message)
                counts_reprocessed <- NA
                proportions_reprocessed <- NA
            })
        } else {
            # Calculate proportions using sweep for raw data
            proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
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

