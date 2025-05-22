<<<<<<< Updated upstream
parse_2023_pereira_nature_nervous <- function(raw = FALSE, align = FALSE) {
    required_pkgs <- c("tidyverse", "readxl", "stringr")
    missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
    if (length(missing_pkgs) > 0) {
        stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
            ". Please install them before running this function.")
    }
    if (!is.logical(raw) || length(raw) != 1) {
        stop("`raw` must be a single logical value (TRUE or FALSE)")
    }
    if (!is.logical(align) || length(align) != 1) {
        stop("`align` must be a single logical value (TRUE or FALSE)")
    }
    
    library(stringr)
    library(tidyverse)

    # ---------- MANAN PROCESSED BELOW ----------

    # localPath <- file.path("/2023_pereira_nature_nervous/")

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

    # ----- Local base directory -----
    local <- file.path("2023_pereira_nature_nervous")

    # ----- File paths -----
    repro_counts_rds_zip <- file.path(local, "PRJNA1033532_dada2_counts.rds.zip")
    repro_tax_zip        <- file.path(local, "PRJNA1033532_dada2_taxa.rds.zip")
    scale_16s_zip        <- file.path(local, "Pereira2023_scale.csv.zip")
    counts_16s_zip       <- file.path(local, "Pereira_2023_16S.csv.zip")
    metadata_16s_zip     <- file.path(local, "Pereira_2023_metadata.csv.zip")
    sra_zip              <- file.path(local, "SraRunTable (38).csv.zip")

    # ----- Initialize everything as NA -----
    counts_original <- NA
    proportions_original <- NA
    tax_original <- NA
    counts_reprocessed <- NA
    proportions_reprocessed <- NA
    counts_reprocessed2 <- NA
    proportions_reprocessed2 <- NA
    tax_reprocessed <- NA
    tax_reprocessed2 <- NA
    scale <- NA
    metadata <- NA

    # ---- scale and metadata -----
    scale         <- read_zipped_table(scale_16s_zip, row.names = NULL) %>% rename(Sample = `Sequencing sample ID`)
    metadata      <- read_zipped_table(metadata_16s_zip, row.names = NULL) %>% rename(Sample = `Sequencing sample ID`)
    sra          <- read_zipped_table(sra_zip, row.names = NULL) %>% rename(Accession = Run, Sample = `Library Name`)
    
    # Merge metadata components
    metadata      <- full_join(scale, metadata, by = "Sample")
    metadata      <- full_join(sra, metadata, by = "Sample")
    
    # Update scale with log transformations
    scale         <- scale %>% 
                    dplyr::select(c("Sample", "Cells/mL")) %>% 
                    mutate(`Cells/mL` = as.numeric(`Cells/mL`)) %>%
                    mutate(log2_FC_cells_ml = ifelse(10^`Cells/mL` > 0, log2(10^`Cells/mL`),NA)) %>% 
                    rename(log10_FC_cells_ml = `Cells/mL`)

    # ------ original counts ------
    counts_original <- read_zipped_table(counts_16s_zip)

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
            original_names <- colnames(counts_original)
            counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
        }

        # ------ proportions from counts ------
        proportions_original <- sweep(counts_original, MARGIN = 1, STATS = rowSums(counts_original), FUN = "/")

    } else {
        proportions_original <- NA
        tax_original <- NA
    }

    if (all(file.exists(repro_counts_rds_zip), file.exists(repro_tax_zip))) {
        # # ----- Reprocessed counts from RDS ZIP -----
        temp_dir <- tempfile("repro")
        dir.create(temp_dir)
        unzipped = unzip(repro_counts_rds_zip, exdir = temp_dir, overwrite = TRUE)
        counts_file <- unzipped[grep("_counts\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(counts_file)) stop("No *_counts.rds file found after unzip")
        counts_reprocessed <- as.data.frame(readRDS(counts_file))

        # ----- rdp16 -----
        if (!file.exists(file.path(local,"rdp16classified.csv.zip"))) {
          if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
              required_pkgs <- c("dada2", "Biostrings")
              missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
              if (length(missing_pkgs) > 0) {
                stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                    ". Please install them before running this function.")
              }
              seqs <- Biostrings::DNAStringSet(colnames(counts_reprocessed))
              rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
              tax_reprocessed2 = make_taxa_label(rdpclassified) 
              write.csv(tax_reprocessed2, file = file.path(local, "rdp16classified.csv"), row.names = TRUE)
            } else {
              stop("RDP 16 file not detected. please install the helperdata/rdp_train_set_16.fa.gz file")
          }
          
          } else {
            tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
        }

        # ----- Taxonomy reprocessed -----
        unzipped = unzip(repro_tax_zip, exdir = temp_dir, overwrite = TRUE)
        tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))
        
        # ----- Convert sequences to lowest rank taxonomy found and update key -----
        tax_reprocessed = make_taxa_label(tax_reprocessed)

        # ----- Convert accessions to sample IDs / Sequences to Taxa -----
        if (!raw) {
            aligned = rename_and_align(counts_reprocessed= counts_reprocessed, metadata=metadata, scale=scale, by_col="Sample", align = align, study_name=basename(local))
            counts_reprocessed = aligned$reprocessed
            counts_reprocessed2 = aligned$reprocessed
            matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
            matched_taxa2 <- tax_reprocessed2$Taxa[match(colnames(counts_reprocessed2), rownames(tax_reprocessed2))]
            colnames(counts_reprocessed) <- matched_taxa
            colnames(counts_reprocessed2) <- matched_taxa2
            counts_reprocessed <- collapse_duplicate_columns_exact(counts_reprocessed)
            counts_reprocessed2 <- collapse_duplicate_columns_exact(counts_reprocessed2)
            original_names <- colnames(counts_reprocessed)
            original_names2 <- colnames(counts_reprocessed2)
            counts_reprocessed <- as.data.frame(lapply(counts_reprocessed, as.numeric), row.names = rownames(counts_reprocessed), col.names = original_names, check.names = FALSE)
            counts_reprocessed2 <- as.data.frame(lapply(counts_reprocessed2, as.numeric), row.names = rownames(counts_reprocessed2), col.names = original_names2, check.names = FALSE)
            proportions_reprocessed2 <- sweep(counts_reprocessed2, 1, rowSums(counts_reprocessed2), '/')
        }

        # proportions reprocessed
        proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
        cleanup_tempfiles(temp_dir)
    }

    if (!raw) {
        counts_original = fill_na_zero_numeric(counts_original)
        proportions_original = fill_na_zero_numeric(proportions_original)
        counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
        proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
        counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
        proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    }


    # ----- Return structured list -----
    return(list(
        counts = list(
            original = counts_original,
            reprocessed = list(rdp19 = counts_reprocessed, rdp16 = counts_reprocessed2)
        ),
        proportions = list(
            original = proportions_original,
            reprocessed = list(rdp19 = proportions_reprocessed, rdp16 = proportions_reprocessed2)
        ),
        tax = list(
            original = tax_original,
            reprocessed = list(rdp19 = tax_reprocessed, rdp16 = tax_reprocessed2)
        ),
        scale = scale,
        metadata = metadata
    ))
=======
parse_2023_pereira_nature_nervous <- function(paths = NULL) {
  required_pkgs <- c("tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  
  require(stringr)
  library(tidyverse)

  localPath <- file.path("/2023_pereira_nature_nervous/")

  # ---------- MANAN PROCESSED BELOW ----------

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

  repro_counts_rds_zip<- paste0(localPath, "PRJNA1033532_dada2_merged_nochim.rds.zip")
  repro_tax_zip       <- paste0(localPath, "PRJNA1033532_dada2_taxonomy_merged.rds.zip")
  scale_16s_zip     <- paste0(local, "Pereira2023_scale.csv.zip")
  counts_16s_zip    <- paste0(local, "Pereira_2023_16S.csv.zip")
  metadata_16s_zip  <- paste0(local, "Pereira_2023_metadata.csv.zip") 

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
>>>>>>> Stashed changes
}