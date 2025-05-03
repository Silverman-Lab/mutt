parse_2024_sinha_natureBiotec_spikeinhumanfecalreproducibility <- function(paths = NULL) {
  required_pkgs <- c("stringr", "tidyverse")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
      stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }

    library(stringr)
    library(tidyverse)

    # ----- Local base directory -----
    local <- file.path("2017_sinha_natureBiotec_spikeinhumanfecalreproducibility")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "SRP047083_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "SRP047083_dada2_taxa.rds.zip")
    counts_16s_zip       <- file.path(local, "mbqc_integrated_otus.tsv.zip")
    metadata_SRA_zip     <- file.path(local, "SraRunTable (38).csv.zip")

    read_zipped_table <- function(zip_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE) {
      if (file.exists(zip_path)) {
        inner_file <- unzip(zip_path, list = TRUE)$Name[1]
        con <- unz(zip_path, inner_file)
        read.table(con, sep = sep, header = header, row.names = row.names, check.names = check.names, stringsAsFactors = FALSE)
      } else {
        warning(paste("File not found:", zip_path))
        return(NA)
      }
    }

    dat <- read_zipped_table(counts_16s_zip, sep="/t")
    metadata <- dat[1:71,]
    metadata <- t(metadata)
    colnames(metadata) <- metadata[1,]
    metadata <- as.data.frame(metadata[-1,])
    
    # dat has metadata (row1-71) and otu counts (row72~) stitched together
    counts <- t(dat[72:nrow(dat),])
    colnames(counts) <- counts[1,]
    counts <- counts[-1,] 
    
    scale <- as.matrix(metadata$nt_count, ncol = 1)
    rownames(scale) <- rownames(metadata)
    
    # Extract taxonomy and the numerical IDs from count table
    tax_list <- strsplit(colnames(counts), "\\|")
    tax <- do.call(rbind, lapply(tax_list, function(taxa) {
      kingdom <- sub("k__", "", taxa[1])
      phylum <- sub("p__", "", taxa[2])
      class <- sub("c__", "", taxa[3])
      order <- sub("o__", "", taxa[4])
      family <- sub("f__", "", taxa[5])
      genus <- sub("g__", "", taxa[6])
      species <- sub("s__", "", taxa[7])
      id <- sub("^.*\\|", "", taxa[length(taxa)])
      c(kingdom, phylum, class, order, family, genus, species, id)
    }))

    tax <- as.data.frame(tax, stringsAsFactors = FALSE)
    colnames(tax) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ID")
    rownames(tax) <- tax$ID
    tax <- tax[, -ncol(tax)]  
    colnames(counts) <- rownames(tax)

    sra = read_zipped_table(metadata_SRA_zip) %>%
      rename(Accession = Run,
             Sample_name = Library.Name
            ) %>%
      mutate(Sample_name = )

    metadata = metadata %>% full_join(sra, by = c(,"Sample_name"))

    
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
    make_taxa_label <- function(df) {
        tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
        prefixes  <- c("k", "p", "c", "o", "f", "g")
        if (!all(tax_ranks %in% colnames(df))) {
        stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
        }
        df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
        x[is.na(x) | trimws(x) == ""] <- "unclassified"
        return(x)
        })
        df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
        for (i in length(tax_ranks):1) {
            if (tax_row[i] != "unclassified") {
            return(paste0(prefixes[i], "_", tax_row[i]))
            }
        }
        return("unclassified")  
        })
        return(df)
    }
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    # accessions to sampleIDs is study specific: IF NEED BE

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
            original = proportions
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
