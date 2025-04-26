parse_2024_jin_pnas_semen <- function(paths = NULLL) {
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
  localPath <- file.path("2024_jin_pnas_semen")

  # ---------- file paths ---------------------
  repro_counts_rds_zip<- file.path(localPath, "PRJNA747100_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- file.path(localPath, "PRJNA747100_dada2_taxonomy_merged.rds.zip")
  
  
  metadata <- as.data.frame(read_table(file.path(localPath, "metadata4.txt")))
  rownames(metadata) <- metadata$`sample-ID`
  metadata <- subset(metadata, select = -`sample-ID`)
  metadata <- metadata[!(rownames(metadata) %in% c("S1", "S13", "SN", "SNTC")), ]
  
  counts <- as.data.frame(read_tsv(paste0(localPath, "table.tsv"), skip = 1))
  rownames(counts) <- counts$`#OTU ID`
  counts <- subset(counts, select = -`#OTU ID`)
  counts <- t(counts)
  counts <- counts[!(rownames(counts) %in% c("S1", "SN", "SNTC")), ]
  
  
  tax <- as.data.frame(read_tsv(paste0(localPath, "taxonomy.tsv")), skip = 1)
  rownames(tax) <- tax$`Feature ID`
  tax <- subset(tax, select = -`Feature ID`)
  
  scale <- as.data.frame(read_table(paste0(localPath, "CFU_data.txt")), skip = 1)
  rownames(scale) <- scale$`sample-ID`
  scale <- subset(scale, select = -`sample-ID`)
  scale <- t(scale)

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
  
  # ----- Reprocessed counts from RDS ZIP -----
  temp_rds            <- tempfile(fileext = ".rds")
  unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)
  rds_file            <- list.files(dirname(temp_rds), pattern = "\\.rds$", full.names = TRUE)[1]
  seqtab_nochim       <- readRDS(rds_file)
  rpt_mat             <- t(seqtab_nochim)
  counts_reprocessed  <- as.data.frame(rpt_mat)
  counts_reprocessed$Sequence <- rownames(counts_reprocessed)
  counts_reprocessed = counts_reprocessed[, c("Sequence", setdiff(names(counts_reprocessed), "Sequence"))]
  rownames(counts_reprocessed) <- paste0("Taxon_", seq_len(nrow(counts_reprocessed)))
  
  # proportions reprocessed
  proportions_reprocessed = counts_reprocessed
  proportions_reprocessed[-1] <- lapply(
    counts_reprocessed[-1],
    function(col) col / sum(col)
  )
  
  # ----- Taxonomy reprocessed -----
  temp_tax <- tempfile(fileext = ".rds")
  unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
  tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
  taxonomy_matrix <- readRDS(tax_file)
  rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
  tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
  tax_reprocessed = tax_table
  
  
  return(list(scale=scale, 
              metadata=metadata, 
              counts=list(
                original = counts, 
                reprocessed = counts_reprocessed
              ),
              proportions=list(
                original = ,
                reprocessed = proportions_reprocessed
              ),
              tax=list(
                original = tax,
                reprocessed = tax_reprocessed
              )
  )
  )
}