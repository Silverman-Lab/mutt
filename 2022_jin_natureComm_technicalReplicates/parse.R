parse_2024_jin_natureComm_technicalReplicates <- function(paths = NULL) {
  required_pkgs <- c("tidyverse", "readxl", "stringr", "readr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages: ", paste(missing_pkgs, collapse = ", "),
      ". Please install them before running this function."
    )
  }

  library(tidyverse)
  library(readxl)
  library(stringr)
  library(readr)
  
  # -----local path ---------------------------
  local <- file.path("2024_jin_natureComm_technicalReplicates")

  # ---------- file paths ---------------------
  repro_counts_zips <- c(
    file.path(local, "PRJNA639639_dada2_merged_nochim.rds.zip"),
    file.path(local, "PRJNA639647_dada2_taxonomy_merged.rds.zip"),
    file.path(local, "PRJNA780331_dada2_taxonomy_merged.rds.zip"),
    file.path(local, "PRJNA780361_dada2_taxonomy_merged.rds.zip")
  )

  repro_tax_zips <- c(
    file.path(local, "PRJNA639639_dada2_taxonomy_merged.rds.zip"),
    file.path(local, "PRJNA639647_dada2_taxonomy_merged.rds.zip"),
    file.path(local, "PRJNA780331_dada2_taxonomy_merged.rds.zip"),
    file.path(local, "PRJNA780361_dada2_taxonomy_merged.rds.zip")
  )
  supp8_zip           <- file.path(local, "Supplementary Data 8.xlsx.zip")


  unzipped_files <- unzip(supp8_zip, exdir = tempdir())
  supp8_xlsx <- unzipped_files[1] 

  dat <- readxl::read_xlsx(supp8_zip, sheet = "Sequencing-determined counts", skip = 1)
  counts <- as.matrix(dat[,2:56])
  rownames(counts) <- dat$`cOTU-ID`
  tax <- dat[,57:62]
  rownames(tax) <- dat$`cOTU-ID`
  scale <- readxl::read_xlsx(supp8_zip, sheet = "Absolute total abundance", skip = 1)
  
  #extract metadata from sample id
  metadata <- as.data.frame(matrix(NA, nrow = nrow(scale), ncol = 5))
  rownames(metadata) <- scale$Sample
  colnames(metadata) <- c("diet", "subject", "location", "sample_type", "technical_replicate")
  
  #first two letter indicate diet type (VD, VS, CE)
  metadata$diet <- substr(scale$Sample, start = 1, stop = 2)
  
  #third letter indicate subject id within group. Note; subject in different group are different mice (i.e., VDa != VSa)
  metadata$subject <- substr(scale$Sample, start = 3, stop = 3)
  
  #sample collected in two different location (distal, proximal)
  metadata$location <- substr(scale$Sample, start = 4, stop = 7)
  
  #sample measured using either filter-residue (cell-sample) and flow-through (ecDNA-sample)
  sampleid <- rownames(metadata)
  metadata$sample_type <- unlist(lapply(strsplit(sampleid, split = "-"), function(x) x[2]))
  metadata$sample_type[is.na(metadata$sample_type)] <- "cell"
  
  #technical replicate if applicable (not every diet group have technical replicates)
  metadata$technical_replicate <- as.numeric(substr(scale$Sample, start = 8, stop = 8))

  # --- Compute proportions from counts ---
  proportions <- counts
  proportions[] <- lapply(
  proportions,
  function(col) {
      if (is.numeric(col)) {
      total <- sum(col, na.rm = TRUE)
      if (total == 0) return(rep(NA, length(col)))
      return(col / total)
      } else {
      return(col)
      }
  }
  )
  
  for (i in seq_along(repro_counts_zips)) {
        # Counts
        temp_rds <- tempfile(fileext = ".rds")
        unzip(repro_counts_zips[i], exdir = dirname(temp_rds), overwrite = TRUE)
        rds_file <- list.files(dirname(temp_rds), pattern = "\\.rds$", full.names = TRUE)[1]
        seqtab_nochim <- readRDS(rds_file)
        rpt_mat <- t(seqtab_nochim)
        counts <- as.data.frame(rpt_mat)
        counts$Sequence <- rownames(counts)
        counts <- counts[, c("Sequence", setdiff(names(counts), "Sequence"))]
        rownames(counts) <- paste0("Taxon_", seq_len(nrow(counts)))
        counts_reprocessed_list[[i]] <- counts

        # Taxonomy
        temp_tax <- tempfile(fileext = ".rds")
        unzip(repro_tax_zips[i], exdir = dirname(temp_tax), overwrite = TRUE)
        tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
        taxonomy_matrix <- readRDS(tax_file)
        rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
        tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
        tax_reprocessed_list[[i]] <- tax_table
    }

  # Combine reprocessed
  counts_reprocessed <- bind_rows(counts_reprocessed_list)
  proportions_reprocessed <- counts_reprocessed
  proportions_reprocessed[-1] <- lapply(
      counts_reprocessed[-1],
      function(col) col / sum(col)
  )
  tax_reprocessed <- bind_rows(tax_reprocessed_list)
  
  return(list(
    scale=scale, 
    metadata=metadata, 
    counts=list(
      original=counts, 
      reprocessed=counts_reprocessed
    ),
    tax=list(
      original=tax,
      reprocessed=tax_reprocessed
    ),
    proportions=list(
      original=proportions,
      reprocessed=proportions_reprocessed
    )
  ))
}


