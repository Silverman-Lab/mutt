parse_2019_morton_naturecommunications_songbird_oral <- function(raw = FALSE, align = FALSE) {
  required_pkgs <- c("stringr", "tidyverse", "matrixStats")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
            stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                ". Please install them before running this function.")
    }
    if (!is.logical(align)) {
      stop("align must be a logical value")
    }
    if (!is.logical(raw)) {
      stop("raw must be a logical value")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2019_morton_naturecommunications_songbird_oral")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "ERP111447_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "ERP111447_dada2_taxa.rds.zip")
    metadata_zip         <- file.path(local, "oral_trimmed_metadata.csv.zip")
    counts_zip           <- file.path(local, "2019_morton_songbird_oral_counts.RDS.zip")
    sra_zip              <- file.path(local, "SraRunTable (40).csv.zip")
    tax_zip              <- file.path(local, "taxonomy.tsv.zip")

    # ----- Initialize -----
    counts = NULL
    proportions = NULL
    tax = NULL
    metadata = NULL
    scale = NULL
    counts_reprocessed = NULL
    proportions_reprocessed = NULL
    tax_reprocessed = NULL

    # ----- Metadata -----
    metadata <- read_zipped_table(metadata_zip, row.names=NULL)
    metadata <- metadata[, !is.na(names(metadata)) & names(metadata) != ""] 

    metadata <- metadata %>%
      mutate(across(c("flow cells/ul 1", "flow cells/ul 2"), as.numeric)) %>%
      mutate(across(c("qPCR cell/ul 1", "qPCR cell/ul 2", "qPCR cell/ul 3"), as.numeric)) %>%
      mutate(across(c("qPCR cell 5 min 1", "qPCR cell 5 min 2", "qPCR cell 5 min 3"), as.numeric)) %>%
      mutate(
        avg_FC_cells_ul = rowMeans(select(., "flow cells/ul 1", "flow cells/ul 2"), na.rm = TRUE),
        sd_FC_cells_ul = apply(select(., "flow cells/ul 1", "flow cells/ul 2"), 1, sd, na.rm = TRUE),
        avg_qpcr_cells_ul = rowMeans(select(., "qPCR cell/ul 1", "qPCR cell/ul 2", "qPCR cell/ul 3"), na.rm = TRUE),
        sd_qpcr_cells_ul = apply(select(., "qPCR cell/ul 1", "qPCR cell/ul 2", "qPCR cell/ul 3"), 1, sd, na.rm = TRUE),
        avg_qpcr_cells_5min = rowMeans(select(., "qPCR cell 5 min 1", "qPCR cell 5 min 2", "qPCR cell 5 min 3"), na.rm = TRUE),
        sd_qpcr_cells_5min = apply(select(., "qPCR cell 5 min 1", "qPCR cell 5 min 2", "qPCR cell 5 min 3"), 1, sd, na.rm = TRUE),
        SampleID = paste0(`participant-timepoint`, ".", Timepoint, ".", treatment)
      ) 

    sra = read_zipped_table(sra_zip, row.names = NULL) %>% rename(Accession = Run)
    sra$SampleID <- paste0(sra$saliva_sample_id, ".", sra$timepoint, ".", sra$processing)
    metadata <- merge(sra, metadata, by = "SampleID")

    scale <- metadata %>% select(SampleID, avg_FC_cells_ul, sd_FC_cells_ul, avg_qpcr_cells_5min, sd_qpcr_cells_5min, avg_qpcr_cells_ul, sd_qpcr_cells_ul) %>% 
      mutate(log2_FC_avg_cells_ul = ifelse(avg_FC_cells_ul > 0, log2(avg_FC_cells_ul), NA)) %>%
      mutate(log10_FC_avg_cells_ul = ifelse(avg_FC_cells_ul > 0, log10(avg_FC_cells_ul), NA)) %>%
      mutate(log2_FC_sd_cells_ul = ifelse(sd_FC_cells_ul > 0, log2(sd_FC_cells_ul), NA)) %>%
      mutate(log10_FC_sd_cells_ul = ifelse(sd_FC_cells_ul > 0, log10(sd_FC_cells_ul), NA)) %>%
      mutate(log2_qpcr_avg_cells_5min = ifelse(avg_qpcr_cells_5min > 0, log2(avg_qpcr_cells_5min), NA)) %>%
      mutate(log10_qpcr_avg_cells_5min = ifelse(avg_qpcr_cells_5min > 0, log10(avg_qpcr_cells_5min), NA)) %>%
      mutate(log2_qpcr_avg_cells_ul = ifelse(avg_qpcr_cells_ul > 0, log2(avg_qpcr_cells_ul), NA)) %>%
      mutate(log10_qpcr_avg_cells_ul = ifelse(avg_qpcr_cells_ul > 0, log10(avg_qpcr_cells_ul), NA)) %>%
      mutate(log2_qpcr_sd_cells_5min = ifelse(sd_qpcr_cells_5min > 0, log2(sd_qpcr_cells_5min), NA)) %>%
      mutate(log10_qpcr_sd_cells_5min = ifelse(sd_qpcr_cells_5min > 0, log10(sd_qpcr_cells_5min), NA)) %>%
      mutate(log2_qpcr_sd_cells_ul = ifelse(sd_qpcr_cells_ul > 0, log2(sd_qpcr_cells_ul), NA)) %>%
      mutate(log10_qpcr_sd_cells_ul = ifelse(sd_qpcr_cells_ul > 0, log10(sd_qpcr_cells_ul), NA))
    

    ## Read Counts
    temp_rds <- tempfile("repro")
    dir.create(temp_rds)
    unzipped = unzip(counts_zip, exdir = dirname(temp_rds), overwrite = TRUE)
    counts_file <- unzipped[grep("2019_morton_songbird_oral_counts\\.RDS$", unzipped, ignore.case = TRUE)][1]
    if (is.na(counts_file)) stop("No 2019_morton_songbird_oral_counts.rds file found after unzip")
    counts <- as.data.frame(readRDS(counts_file)) %>% t() %>% as.data.frame()
    
    ## Taxonomy Information
    raw_tax <- read_zipped_table(tax_zip, sep="\t", row.names = NULL)
    tax <- raw_tax %>%
        mutate(taxonomy = str_replace_all(Taxon, "\\s+", "")) %>%
        mutate(ogtaxonomy = taxonomy) %>%
        separate(
            taxonomy,
            into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
            sep = ";",
            extra = "drop",
            fill = "right"
        ) %>%
        mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                     ~gsub("^[kpcofgs]__", "", .))) %>%
        mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                     ~gsub("\\[|\\]|\\(|\\)", "", .))) %>% 
        mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                     ~gsub("_+", "_", .))) %>% 
        mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus, Species), 
                     ~ifelse(. == "" | . == "__" | is.na(.), "unclassified", .)))
    tax = make_taxa_label(tax)
    tax <- tax[,c("Feature ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Taxa", "ogtaxonomy")]
    row.names(tax) <- tax$`Feature ID`

    if (!raw) {
        aligned = rename_and_align(counts_original = counts, metadata=metadata, scale=scale, by_col = "anonymized_name", align = align, study_name = basename(local))
        counts <- aligned$counts_original
        matched_taxa <- tax$Taxa[match(colnames(counts), rownames(tax))]
        colnames(counts) <- matched_taxa
        counts = collapse_duplicate_columns_exact(counts)
        original_names <- colnames(counts)
        counts <- as.data.frame(lapply(counts, as.numeric), row.names = rownames(counts), col.names = original_names, check.names = FALSE)
    }
    proportions <- sweep(counts, 1, rowSums(counts), FUN = "/")
    
    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
        unzipped = unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(counts_file))

        # ----- Taxonomy reprocessed -----
        unzipped = unzip(repro_tax_zip, exdir = dirname(temp_rds), overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))
        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col = "SampleID", align = align, study_name = basename(local))
            counts_reprocessed <- aligned$reprocessed
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            colnames(counts_reprocessed) <- matched_taxa
            counts_reprocessed = collapse_duplicate_columns_exact(counts_reprocessed)   
            original_names <- colnames(counts_reprocessed)
            counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
        }

        # proportions reprocessed
        proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
    }

    if (!raw) {
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        counts = fill_na_zero_numeric(counts)
        proportions = fill_na_zero_numeric(proportions)
    }   

    cleanup_tempfiles(temp_rds)

    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts,
            reprocessed = counts_reprocessed
        ),
        proportions = list(
            original = proportions,
            reprocessed = proportions_reprocessed
        ),
        tax = list(
            original = tax,
            reprocessed = tax_reprocessed
        ),
        scale = scale,
        metadata = metadata
    ))
}

