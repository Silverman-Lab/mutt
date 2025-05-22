<<<<<<< Updated upstream
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
    nishijima2024counts  <- file.path(local, "Vandeputte_2017_16S.zip")

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

    metadata <- full_join(metadata, sra, by = "Accession") %>% rename(Sample = Sample.x)

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
        mutate(log2_FC_cell_g_frozen = ifelse(`Average cell count (per gram of frozen feces)`>0, log2(`Average cell count (per gram of frozen feces)`), NA)) %>%
        mutate(log10_FC_cell_g_frozen = ifelse(`Average cell count (per gram of frozen feces)`>0, log10(`Average cell count (per gram of frozen feces)`), NA))  %>% 
        mutate(log2_FC_cell_g_fresh = ifelse(`Average cell count (per gram of fresh feces)`>0, log2(`Average cell count (per gram of fresh feces)`), NA)) %>%
        mutate(log10_FC_cell_g_fresh = ifelse(`Average cell count (per gram of fresh feces)`>0, log10(`Average cell count (per gram of fresh feces)`), NA)) %>%
        as.data.frame()  %>%
        select(-Day) %>%
        mutate(cohort = Sample %>% str_extract("^[A-Z]+"))

    # ----- Original counts from CSV.zip -----
    if (file.exists(orig_counts_zip)) {
        orig_csv <- unzip(orig_counts_zip, list = TRUE)$Name[1]
        orig_con <- unz(orig_counts_zip, orig_csv)
        orig_mat <- read.csv(orig_con, row.names = 1, check.names = FALSE)

        if (!file.exists(file.path(local,"otu_taxonomy_rdp.csv.zip"))) {
            required_pkgs <- c("dada2")
            missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
            if (length(missing_pkgs) > 0) {
            stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                ". Please install them before running this function.")
            }
            
            asv_seqs <- colnames(orig_mat)
            names(asv_seqs) <- asv_seqs
            taxa <- dada2::assignTaxonomy(asv_seqs,"helperdata/rdp_19_toGenus_trainset.fa.gz",multithread = TRUE)
            orig_tax_rdp <- as.data.frame(taxa, stringsAsFactors=FALSE)
            orig_tax_rdp$ASV <- rownames(orig_tax_rdp)
            write.csv(orig_tax_rdp, file.path(local,"otu_taxonomy_rdp.csv"), row.names=FALSE)
        } else {
            orig_tax_rdp_zip <- file.path(local, "otu_taxonomy_rdp.csv.zip")
            orig_tax_rdp <- read_zipped_table(orig_tax_rdp_zip, row.names=NULL,check.names = FALSE)
        }

        if (!file.exists(file.path(local,"otu_taxonomy_silva.csv.zip"))) {
            required_pkgs <- c("dada2")
            missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
            if (length(missing_pkgs) > 0) {
            stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                ". Please install them before running this function.")
            }
            asv_seqs <- colnames(orig_mat)
            names(asv_seqs) <- asv_seqs
            if (file.exists("helperdata/silva_nr99_v138.1_train_set.fa.gz")) {
                taxa <- dada2::assignTaxonomy(asv_seqs,"helperdata/silva_nr99_v138.1_train_set.fa.gz",multithread = TRUE)
            } else {
                 stop("Silva classifier file not detected: helperdata/silva_nr99_v138.1_train_set.fa.gz, Please download silva_nr99_v138.1_train_set.fa.gz before running this function.")
            }
            orig_tax_silva <- as.data.frame(taxa, stringsAsFactors=FALSE)
            orig_tax_silva$ASV <- rownames(orig_tax_silva)
            write.csv(orig_tax_silva, file.path(local,"otu_taxonomy_silva.csv"), row.names=FALSE)
        } else {
            orig_tax_silva_zip <- file.path(local, "otu_taxonomy_silva.csv.zip")
            orig_tax_silva <- read_zipped_table(orig_tax_silva_zip, row.names=1,check.names = FALSE) %>% rownames_to_column("ASV")
        }

        if (!file.exists(file.path(local,"otu_taxonomy_rdp16.csv.zip"))) {
            if (file.exists(file.path("helperdata/rdp_train_set_16.fa.gz"))) {
                required_pkgs <- c("dada2", "Biostrings")
                missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
                if (length(missing_pkgs) > 0) {
                    stop("RDP classifier detected. Missing required packages: ", paste(missing_pkgs, collapse = ", "),
                        ". Please install them before running this function.")
                }
                seqs <- Biostrings::DNAStringSet(colnames(orig_mat))
                rdpclassified <- dada2::assignTaxonomy(seqs, file.path("helperdata/rdp_train_set_16.fa.gz"), multithread=TRUE) %>% as.data.frame()
                orig_tax_rdp16 <- as.data.frame(rdpclassified, stringsAsFactors=FALSE)
                orig_tax_rdp16$ASV <- rownames(orig_tax_rdp16)
                write.csv(orig_tax_rdp16, file.path(local,"otu_taxonomy_rdp16.csv"), row.names=FALSE)
            } else {
                stop("RDP classifier file not detected: helperdata/rdp_train_set_16.fa.gz, Please download rdp_train_set_16.fa.gz before running this function.")
            }
        } else {
            orig_tax_rdp16 <- read_zipped_table(file.path(local, "otu_taxonomy_rdp16.csv.zip"), row.names = "ASV")
        }

        tax_original_rdp16 = make_taxa_label(orig_tax_rdp16)


        original_tax_classified <- as.data.frame(merge(orig_tax_rdp, orig_tax_silva, by="ASV", all=TRUE, suffixes = c("_rdp", "_silva")))
        orig_tax <- data.frame(
            Sequence = original_tax_classified$ASV,
            Kingdom  = original_tax_classified$Kingdom_rdp,
            Phylum   = original_tax_classified$Phylum_rdp,
            Class    = original_tax_classified$Class_rdp,
            Order    = original_tax_classified$Order_rdp,
            Family   = original_tax_classified$Family_rdp,
            Genus    = original_tax_classified$Genus_rdp,
            Kingdom_silva = original_tax_classified$Kingdom_silva,
            Phylum_silva  = original_tax_classified$Phylum_silva,
            Class_silva   = original_tax_classified$Class_silva,
            Order_silva   = original_tax_classified$Order_silva,
            Family_silva  = original_tax_classified$Family_silva,
            Genus_silva   = original_tax_classified$Genus_silva
        )

        tax_original <- make_taxa_label(orig_tax)
        tax_original <- tax_original %>% rename(Kingdom_rdp = Kingdom, Phylum_rdp = Phylum, Class_rdp = Class, Order_rdp = Order, Family_rdp = Family, Genus_rdp = Genus)
        rownames(tax_original) <- tax_original$Sequence
        counts_original <- orig_mat


        if (!raw) {
            aligned <- rename_and_align(counts_original = counts_original, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
            counts_original <- aligned$counts_original
            counts_original2 <- aligned$counts_original
            matched_taxa <- tax_original$Taxa[match(colnames(counts_original), rownames(tax_original))]
            matched_taxa2 <- tax_original_rdp16$Taxa[match(colnames(counts_original), rownames(tax_original))]
            colnames(counts_original) <- matched_taxa
            colnames(counts_original2) <- matched_taxa2
            counts_original <- collapse_duplicate_columns_exact(counts_original)
            counts_original2 <- collapse_duplicate_columns_exact(counts_original2)
            original_names <- colnames(counts_original)
            original_names2 <- colnames(counts_original2)
            counts_original <- as.data.frame(lapply(counts_original, as.numeric), row.names = rownames(counts_original), col.names = original_names, check.names = FALSE)
            counts_original2 <- as.data.frame(lapply(counts_original2, as.numeric), row.names = rownames(counts_original2), col.names = original_names2, check.names = FALSE)
            proportions_original2 <- sweep(counts_original2, MARGIN = 1, STATS = rowSums(counts_original2), FUN = "/")
        }

        proportions_original <- sweep(counts_original, MARGIN = 1, STATS = rowSums(counts_original), FUN = "/")

    } else if (!file.exists(orig_counts_zip) && file.exists(orig_prop_zip)) {
        counts_original <- NA
        counts_original2 <- NA
        proportions_original2 <- NA
        tax_original_rdp16 = NA
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
        counts_original2 <- NA
        proportions_original2 <- NA
        tax_original_rdp16 = NA
    }

    # ----- Nishijima 2024 counts from CSV ZIP -----
    if (file.exists(nishijima2024counts)) {
        nishijima2024counts <- read_zipped_table(nishijima2024counts, sep="\t", row.names = 1, check.names = FALSE)
        if (!raw) {
            aligned <- rename_and_align(counts_reprocessed = nishijima2024counts, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
            nishijima2024counts <- aligned$reprocessed
            nishijima2024counts <- as.data.frame(lapply(nishijima2024counts, as.numeric), row.names = rownames(nishijima2024counts), col.names = colnames(nishijima2024counts), check.names = FALSE)
        }
        proportions_nishijima2024 <- sweep(nishijima2024counts, 1, rowSums(nishijima2024counts), '/')
        tax_nishijima2024 <- data.frame(Taxa = colnames(nishijima2024counts), stringsAsFactors = FALSE)
    }

    # ----- Reprocessed counts from RDS ZIP -----
    if (all(file.exists(repro_counts_rds_zip))) {
        temp_dir <- tempfile("repro")   
        dir.create(temp_dir)
        
        # Unzip counts file
        unzipped =unzip(repro_counts_rds_zip, exdir = temp_dir)
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
            NULL
          }
          
          } else {
            tax_reprocessed2 <- read_zipped_table(file.path(local, "rdp16classified.csv.zip"), row.names = 1)
        }
        
        unzipped = unzip(repro_tax_zip, exdir = temp_dir)
        tax_file <- unzipped[grep("_taxa\\.rds$",   unzipped, ignore.case = TRUE)][1]
        if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
        tax_reprocessed <- as.data.frame(readRDS(tax_file))
        
        # Convert sequences to lowest rank taxonomy found and update key
        tax_reprocessed = make_taxa_label(tax_reprocessed)
        
        # Convert accessions to sample IDs / Sequences to Taxa
        if (!raw) {
            tryCatch({
                aligned <- rename_and_align(counts_reprocessed = counts_reprocessed, metadata = metadata, scale = scale, by_col = "Sample", align = align, study_name = basename(local))
                counts_reprocessed <- aligned$reprocessed
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
                proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
            }, error = function(e) {
                warning("Error in processing reprocessed data: ", e$message)
                counts_reprocessed <- NA
                proportions_reprocessed <- NA
            })
        } else {
            # Calculate proportions using sweep for raw data
            proportions_reprocessed <- sweep(counts_reprocessed, 1, rowSums(counts_reprocessed), '/')
            proportions_reprocessed2 = NA
        }
        cleanup_tempfiles(temp_dir)
    } else {
        counts_reprocessed <- NA
        tax_reprocessed <- NA
        proportions_reprocessed <- NA
        counts_reprocessed2 <- NA
        tax_reprocessed2 <- NA
        proportions_reprocessed2 <- NA
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
            if (!is.na(proportions_reprocessed2)[1]) {
                proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
            }
            if (!is.na(counts_reprocessed2)[1]) {
                counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
            }
        }, error = function(e) {
            warning("Error filling NA values: ", e$message)
        })
    }

    # ----- Return all -----
    return(list(
      counts      = list(
        original    = list(
          rdp19 = counts_original,
          rdp16 = counts_original2
        ),
        nishijima2024 = nishijima2024counts,
        reprocessed = list(
          rdp19 = counts_reprocessed,
          rdp16 = counts_reprocessed2
        )
      ),
      proportions = list(
        original    = list(
          rdp19 = proportions_original,
          rdp16 = proportions_original2
        ),
        nishijima2024 = proportions_nishijima2024,
        reprocessed = list(
          rdp19 = proportions_reprocessed,
          rdp16 = proportions_reprocessed2
        )
      ),
      tax         = list(
        original = list(
          rdp19 = tax_original,
          rdp16 = tax_original_rdp16
        ),
        nishijima2024 = tax_nishijima2024,
        reprocessed = list(
          rdp19 = tax_reprocessed,
          rdp16 = tax_reprocessed2
        )
      ),
      scale       = scale,
      metadata    = metadata
    ))
}
=======
parse_2017_vandeputte_nature_flow <- function() {
  required_pkgs <- c("tibble", "tidyverse", "readxl")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
  library(tibble)
  library(tidyverse)
  library(readxl)

  # ----- File paths -----
  metadata_zip        <- "2017_vandeputte_nature_flow/Vandeputte_2017_metadata.csv.zip"
  metadata_two_zip    <- "2017_vandeputte_nature_flow/cellcountstotal.csv.zip"
  orig_counts_zip     <- "2017_vandeputte_nature_flow/OTU_nochim.zip"
  orig_tax_rdp_zip    <- "2017_vandeputte_nature_flow/otu_taxonomy_rdp.csv.zip"
  orig_tax_silva_zip  <- "2017_vandeputte_nature_flow/otu_taxonomy_silva.csv.zip"
  orig_prop_zip       <- NA
  repro_counts_rds_zip<- "2017_vandeputte_nature_flow/PRJEB21504_dada2_merged_nochim.rds.zip"
  repro_tax_zip       <- "2017_vandeputte_nature_flow/PRJEB21504_dada2_taxonomy_merged.rds.zip"

  # ----- Metadata and Scale -----
  # Read scale
  count_csv    <- unzip(metadata_two_zip, list = TRUE)$Name[1]
  count_con    <- unz(metadata_two_zip, count_csv)
  metadata_two <- read.csv(count_con) %>% as.data.frame()

  # Read metadata
  meta_csv     <- unzip(metadata_zip, list = TRUE)$Name[1]
  meta_con     <- unz(metadata_zip, meta_csv)
  metadata_df  <- read.csv(meta_con) %>% as.data.frame()

  # Join and select
  df <- left_join(
    metadata_df,
    metadata_two %>% select(Sample, SampleID) %>% rename(Accession = SampleID),
    by = c("Individual" = "Sample")
  )
  metadata = df %>% select(
    Individual, Cohort, Day,
    `Health status`, accession = Accession, Enterotype
  )
  scale = df %>% select(
    Individual, Day, accession = Accession,
    `Average cell count (per gram of fresh feces)`,
    `STDEV cell count (per gram of fresh feces)`,
    `Average cell count (per gram of frozen feces)`,
    `STDEV cell count (per gram of frozen feces)`
  )

  # ----- Original counts from CSV.zip -----
  if (file.exists(orig_counts_zip)) {
    orig_csv <- unzip(orig_counts_zip, list = TRUE)$Name[1]
    orig_con <- unz(orig_counts_zip, orig_csv)
    orig_mat <- read.csv(orig_con, row.names = 1, check.names = FALSE)
    counts_original <- as.data.frame(orig_mat)
    counts_original$Sequence <- rownames(counts_original)
    counts_original <- counts_original[, c("Sequence", setdiff(names(counts_original), "Sequence"))]
    rownames(counts_original) <- paste0("Taxon_", seq_len(nrow(counts_original)))
  } else {
    counts_original <- NA
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

  # --- Original taxonomies ---
  read_taxonomy_zip <- function(zip_path) {
    if (!file.exists(zip_path)) return(NA)
    zip_contents <- unzip(zip_path, list = TRUE)
    csv_name <- zip_contents$Name[1]
    con <- unz(zip_path, csv_name)
    tax_df <- read.csv(con, row.names = 1, check.names = FALSE)
    tax_df$Sequence <- rownames(tax_df)
    tax_df <- tax_df[, c("Sequence", setdiff(names(tax_df), "Sequence"))]
    rownames(tax_df) <- paste0("Taxon_", seq_len(nrow(tax_df)))
    as_tibble(tax_df, rownames = "Taxon")
  }

  tax_original_rdp   <- read_taxonomy_zip(orig_tax_rdp_zip)
  tax_original_silva <- read_taxonomy_zip(orig_tax_silva_zip)

  # ----- Taxonomy reprocessed -----
  temp_tax <- tempfile(fileext = ".rds")
  unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)
  tax_file <- list.files(dirname(temp_tax), pattern = "\\.rds$", full.names = TRUE)[1]
  taxonomy_matrix <- readRDS(tax_file)
  rownames(taxonomy_matrix) <- paste0("Taxon_", seq_len(nrow(taxonomy_matrix)))
  tax_table <- as_tibble(taxonomy_matrix, rownames = "Taxon")
  tax_reprocessed <- tax_table

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
      original_rdp    = tax_original_rdp,
      original_silva  = tax_original_silva,
      reprocessed = tax_reprocessed
    ),
    scale       = scale,
    metadata    = metadata
  ))
}
>>>>>>> Stashed changes
