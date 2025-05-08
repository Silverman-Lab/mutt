parse_2017_sinha_natureBiotec_spikeinhumanfecalreproducibility <- function(raw = FALSE, align = FALSE) {
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
  local <- file.path("2017_sinha_natureBiotec_spikeinhumanfecalreproducibility")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "SRP047083_dada2_counts.rds.zip")
  repro_tax_zip        <- file.path(local, "SRP047083_dada2_taxa.rds.zip")
  counts_16s_zip       <- file.path(local, "mbqc_integrated_otus.tsv.zip")
  metadata_SRA_zip     <- file.path(local, "SraRunTable (38).csv.zip")
  metadatareadin       <- file.path(local, "metadata.csv.zip")
  countsreadin         <- file.path(local, "counts.csv.zip")
  proportionsreadin    <- file.path(local, "proportions.csv.zip")
  scalereadin          <- file.path(local, "scale.csv.zip")
  taxreadin            <- file.path(local, "tax.csv.zip")

  # ------ original counts, scale, metadata, proportions, taxa -------

  ## DONT DELETE --- THIS IS THE PROCESSING CODE, BUT I SAVED THE RESULT FOR EFFICIENCY -- huge files

  # dat <- read_zipped_table(counts_16s_zip, sep = "\t")
  # metadata <- as.data.frame(t(dat[   1:71, , drop = FALSE]),
  #                           stringsAsFactors = FALSE)
  # counts   <- as.data.frame(t(dat[-(1:71), , drop = FALSE]),
  #                           stringsAsFactors = FALSE)
  # rm(dat)
  # gc()

  # scale <- metadata %>%
  #           pull(nt_count) %>%             
  #           as.numeric() %>%                  
  #           matrix(ncol = 1) %>%              
  #           `rownames<-`(rownames(metadata))  

  # tax <- tibble(raw = colnames(counts)) %>%
  #   separate(raw,
  #           into  = c("Kingdom","Phylum","Class","Order",
  #                     "Family","Genus","Species","ID"),
  #           sep   = "\\|",
  #           fill  = "right",
  #           extra = "merge") %>%
  #   mutate(across(Kingdom:Species, ~ sub("^[a-z]__", "", .))) %>%
  #   make_taxa_label() %>%
  #   as.data.frame(stringsAsFactors = FALSE)

  # counts <- as.data.frame(counts, stringsAsFactors = FALSE)
  # if ("Sample" %in% colnames(counts)) {
  #   rownames(counts) <- counts$Sample
  #   counts$Sample     <- NULL
  # }
  # colnames(counts) <- tax$Taxa
  # counts <- t(
  #   rowsum(
  #     t(counts),
  #     group = colnames(counts),
  #     reorder = FALSE  
  #   )
  # )
  # counts <- as.data.frame(counts)
  # tax <- tax[colnames(counts), ]
  # proportions <- sweep(counts, MARGIN = 1,STATS  = rowSums(counts), FUN = "/")

  # # Save counts
  # write.csv(counts, file = file.path(local, "counts.csv"), row.names = TRUE)

  # # Save proportions
  # write.csv(proportions, file = file.path(local, "proportions.csv"), row.names = TRUE)

  # # Save scale (1-column matrix)
  # scale_df <- as.data.frame(scale)
  # write.csv(scale_df, file = file.path(local, "scale.csv"), row.names = TRUE)

  # # Save taxonomy table
  # write.csv(tax, file = file.path(local, "tax.csv"), row.names = TRUE)

  ## DONT DELETE --- PROCESSING CODE, BUT SAVED FOR EFFICIENCY

  # ----- Load in processed sheets from code above -----
  metadata     <- read_zipped_table(metadatareadin, row.names = NULL)
  counts       <- read_zipped_table(countsreadin)
  proportions  <- read_zipped_table(proportionsreadin)
  scale        <- read_zipped_table(scalereadin)
  tax          <- read_zipped_table(taxreadin)
  sra          <- read_zipped_table(metadata_SRA_zip, row.names = NULL)

  sra = read_zipped_table(metadata_SRA_zip, row.names = NULL) 
  sra = sra %>%
    rename(
      Accession = Run,
      Sample_name = `Library Name`
      ) 
  metadata <- tryCatch(
    {
      metadata %>% full_join(sra, by = "Sample_name")
    },
    error = function(e) {
      warning("Joining metadata and sra failed: ", e$message)
      metadata  # return the original if there was an error
    }
  )

  if (!raw) {
    align <- rename_and_align(counts_original = counts, proportions_original = proportions, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
    counts = align$counts_original
    proportions = align$proportions_original
  }

  # ----- Calculate scale factor for spike-ins -----
  spike_taxa <- c() # I need to figure out which spike-in taxa were used per sample.

  if (length(spike_taxa) > 0) {
    spikein_counts <- counts[, colnames(counts) %in% spike_taxa, drop = FALSE]
    observed_spike_in_reads <- rowSums(spikein_counts)
    scale_factor = scale / observed_spike_in_reads
  } else {
    scale_factor = scale
  }
  
  if (all(file.exists(repro_counts_rds_zip))) {
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
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
      align <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
      counts_reprocessed = align$reprocessed
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
  }

  if (length(spike_taxa) > 0) {
    spikein_counts <- counts_reprocessed[, colnames(counts_reprocessed) %in% spike_taxa, drop = FALSE]
    observed_spike_in_reads <- rowSums(spikein_counts)
    reprocessed_scale_factor = scale / observed_spike_in_reads
  } else {
    reprocessed_scale_factor = scale
  }

  if (!raw) {
    counts = fill_na_zero_numeric(counts)
    counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
    proportions = fill_na_zero_numeric(proportions)
    proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
  }

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
      scale = list(
        recordedspikein        = scale,
        originalscalefactor    = scale_factor,
        reprocessedscalefactor = reprocessed_scale_factor
      ),
      metadata = metadata
  ))
}
