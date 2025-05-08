parse_2019_morton_naturecommunications_songbird_oral <- function(raw = FALSE, align = align) {

  # THIS STUDY IS THE SAME AS MAROTZ ET AL. AND CAN BE ARCHIVED WITHOUT FURTHER WORK.
  #
  #
  #
  # IF YOU SEE THIS, PLEASE FIGURE OUT BEST WAY TO ARCHIVE THIS> MAYBE COMBINE THE DATASHEETS AND RELABEL WITH MORTON IN THE
  # MAROTZ DIRECTORY? I DONT KNOW.
  #
  #
  #
  #
  required_pkgs <- c("stringr", "tidyverse")
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
    repro_counts_rds_zip <- file.path(local, "PRJNA1233249_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA1233249_dada2_taxa.rds.zip")
    #scale_zip            <- file.path(local, "PMID-28743816_samples-v1.csv.zip")
    metadata_zip         <- file.path(local, "oral_trimmed_metadata.csv")
    counts_zip           <- file.path(local, "2019_morton_songbird_oral_counts.RDS.zip")

    ## Get metadata, add average flowcount, round to avoid fraction of a cell
    metadata <- read_zipped_table(metadata_zip)
    metadata$avg_norm_flowcounts <- rowMeans(metadata[,c("flow.cells.ul.1", "flow.cells.ul.2")])
    metadata$sd_norm_flowcounts <- rowSds(as.matrix(metadata[,c("flow.cells.ul.1", "flow.cells.ul.2")]))
    scale <- scale %>% select(ID, avg_norm_flowcounts, sd_norm_flowcounts) %>% 
                        mutate(log2_flowcounts = ifelse(avg_norm_flowcounts > 0, log2(avg_norm_flowcounts), NA)) %>%
                        mutate(log10_flowcounts = ifelse(avg_norm_flowcounts > 0, log10(avg_norm_flowcounts), NA)) %>%
                        mutate(log2_flowcounts_sd = ifelse(sd_norm_flowcounts > 0, log2(sd_norm_flowcounts), NA)) %>%
                        mutate(log10_flowcounts_sd = ifelse(sd_norm_flowcounts > 0, log10(sd_norm_flowcounts), NA))
    names(scale) <- row.names(metadata)

    ## Read Counts
    temp_rds <- tempfile(fileext = ".rds")
    unzip(counts_zip, exdir = dirname(temp_rds), overwrite = TRUE)

    rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
    if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
    counts<- as.data.frame(readRDS(rds_files[1]))
    if (!raw) {
        align = rename_and_align(counts_original= counts, metadata, scale, by_col = "ID", align = align, study_name = basename(local))
        counts <- align$counts_original
    }
    proportions <- apply(counts, 2, function(col) col/sum(col))

    ## Taxonomy Information
    raw_tax <- read.csv("taxonomy.tsv", sep="\t")
    tax <- raw_tax %>%
        mutate(Taxon=str_trim(Taxon)) %>%
        separate(
        Taxon,
        into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
        sep = "\\s*;\\s*",
        extra = "drop",
        fill = "right"
    )
    row.names(tax) <- tax$Feature.ID
    tax <- tax[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]

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
        align = rename_and_align(counts_reprocessed = counts_reprocessed, metadata, scale, by_col = "ID", align = align, study_name = basename(local))
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

