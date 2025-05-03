parse_2023_pereira_nature_nervous <- function(paths = NULL) {
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  
  require(stringr)
  library(tidyverse)

<<<<<<< Updated upstream
<<<<<<< Updated upstream
  # ---------- MANAN PROCESSED BELOW ----------

  # localPath <- "/2023_pereira_nature_nervous/"

=======
=======
>>>>>>> Stashed changes
  localPath <- file.path("/2023_pereira_nature_nervous/")

  # ---------- MANAN PROCESSED BELOW ----------

<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes
  # decompressed_file <- gunzip(paste0(localPath, "all-relevant-data.xlsx.gz"), remove = FALSE)
  
  # totalCounts <- readxl::read_xlsx(decompressed_file, sheet = "Sup. Table 1", skip = 1)
  # totalCounts <- totalCounts[, c("time (hours)", "condition", "medium", "Cells/mL")]
  # totalCounts$index <- paste(totalCounts[["time (hours)"]], totalCounts$condition, totalCounts$medium, sep = "_")
  # totalCounts$replicate <- ave(totalCounts$index, totalCounts$index, FUN = seq_along)
  # totalCounts$index <- paste(totalCounts$index, totalCounts$replicate, sep = "_")
  # totalCounts <- totalCounts[, c("index", "Cells/mL")]
  
  # sampleMetadata <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 2", skip = 3)
  
  # # Step 1: Process sampleMetadata
  # sampleMetadata <- sampleMetadata %>%
  #   mutate(index = paste(`time (hours)`, condition, medium, sep = "_")) %>%
  #   group_by(index) %>%
  #   mutate(replicate = seq_along(index)) %>%
  #   ungroup() %>%
  #   mutate(index = paste(index, replicate, sep = "_"))
  
  # # Step 2: Intersect with totalCounts
  # scale_data <- sampleMetadata %>%
  #   inner_join(totalCounts, by = "index") %>%
  #   dplyr::select("Sequencing sample ID", `Cells/mL`)
  
  # # Step 3: Create scale matrix
  # scale <- scale_data %>%
  #   pivot_wider(names_from = "Sequencing sample ID", values_from = `Cells/mL`) %>%
  #   as.matrix()
  
  # relativeAbundances <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 3", skip = 3)
  # absoluteAbundances <- readxl::read_xlsx(paste0(localPath, "all-relevant-data.xlsx"), sheet = "Sup. Table 4", skip = 3)
  
  # absoluteAbundances <- absoluteAbundances %>%
  #   rename(ASV_ID = `...1`)
  
  # relativeAbundances <- relativeAbundances %>%
  #   rename(ASV_ID = `...1`)
  
  # countsDataFile <- unzip(paste0(localPath, "DADA2_counts_as_matrix_Drugs.txt.zip"))
  
  # countsData <- read.table(countsDataFile, header = TRUE, sep = "\t")
  
  # counts <- countsData %>% column_to_rownames("X") %>% as.matrix()
  
  # proportions <- relativeAbundances %>%
  #   dplyr::select(-c(Phylum, Class, Order, Family, Genus)) %>%
  #   column_to_rownames("ASV_ID") %>%
  #   as.matrix()
  
  # metadata <- sampleMetadata %>%
  #   column_to_rownames("Sequencing sample ID")
  
  # tax <- absoluteAbundances %>%
  #   dplyr::select(ASV_ID, Phylum, Class, Order, Family, Genus) %>%
  #   column_to_rownames("ASV_ID")
  
  # unlink(decompressed_file)
  # unlink(countsDataFile)

  # ------ MAXWELL PROCESSED BELOW ------------

<<<<<<< Updated upstream
<<<<<<< Updated upstream
  # ----- Local base directory -----
  local <- file.path("2023_pereira_nature_nervous")

  # ----- File paths -----
  repro_counts_rds_zip <- file.path(local, "PRJNA1033532_dada2_merged_nochim.rds.zip")
  repro_tax_zip        <- file.path(local, "PRJNA1033532_dada2_taxonomy_merged.rds.zip")
  scale_16s_zip        <- file.path(local, "Pereira2023_scale.csv.zip")
  counts_16s_zip       <- file.path(local, "Pereira_2023_16S.csv.zip")
  metadata_16s_zip     <- file.path(local, "Pereira_2023_metadata.csv.zip")

=======
=======
>>>>>>> Stashed changes
  repro_counts_rds_zip<- paste0(localPath, "PRJNA1033532_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- paste0(localPath, "PRJNA1033532_dada2_taxonomy_merged.rds.zip")
  scale_16s_zip     <- paste0(local, "Pereira2023_scale.csv.zip")
  counts_16s_zip    <- paste0(local, "Pereira_2023_16S.csv.zip")
  metadata_16s_zip  <- paste0(local, "Pereira_2023_metadata.csv.zip") 
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

  read_zipped_csv <- function(zip_path) {
    if (file.exists(zip_path)) {
        csv_file <- unzip(zip_path, list = TRUE)$Name[1]
        read.csv(unz(zip_path, csv_file), row.names = 1, check.names = FALSE)
    } else {
        warning(paste("File not found:", zip_path))
        return(NA)
    }
  }

  # ----- Initialize everything as NA -----
  counts_original <- NA
  proportions_original <- NA
  tax_original <- NA
  counts_reprocessed <- NA
  proportions_reprocessed <- NA
  tax_reprocessed <- NA

  # ------ original counts ------
  counts_original <- read_zipped_csv(counts_16s_zip)

  if (!is.na(counts_original)[1]) {
    original_taxa <- colnames(counts_original)
    taxon_ids <- paste0("Taxon_", seq_len(nrow(counts_original)))
    colnames(counts_original) <- taxon_ids

    # Create taxa mapping data frame
    tax_original <- data.frame(
      Taxon = taxon_ids,
      Original_Taxa = original_taxa,
      stringsAsFactors = FALSE
    )

    # ------ proportions from counts ------
    proportions_original <- t(counts_original)
    proportions_original[] <- lapply(
      proportions_original,
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
  } else {
    proportions_original <- NA
    tax_original <- NA
  }

  # ---- scale and metadata -----
  scale     <- read_zipped_csv(scale_16s_zip)
  metadata  <- read_zipped_csv(metadata_16s_zip)

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
  
  
  return(list(
    counts = list(
      original = counts_original, 
      reprocessed = counts_reprocessed
    )
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