parse_2019_vieirasilva_naturemicrobiology_pscibd <- function(raw = FALSE, align = FALSE, entrez_key = NULL) {
  # if (is.null(entrez_key) || !nzchar(entrez_key)) {
  #   message("No NCBI Entrez API key provided.")
  #   response <- readline(prompt = "Would you like to enter an NCBI API key? (y/n): ")
  #   if (tolower(substr(response, 1, 1)) == "y") {
  #     entrez_key <- readline(prompt = "Please enter your NCBI API key: ")
  #     if (!nzchar(entrez_key)) {
  #       message("No key entered. Proceeding without API key.")
  #     }
  #   } else {
  #     message("Proceeding without API key.")
  #     message(
  #     "You should provide a valid NCBI Entrez API key via the `entrez_key` argument.\n",
  #     "To create or retrieve your key:\n",
  #     "  1) Visit: https://www.ncbi.nlm.nih.gov/account/\n",
  #     "  2) Sign in or create an account\n",
  #     "  3) Go to 'Settings' → 'API Keys' and generate a new key\n",
  #     "  4) Provide that key to this function, e.g.:\n",
  #     "     parse_2019_viera-silva_naturemicrobiology_pscibd(entrez_key = \"YOUR_KEY_HERE\")"
  #     )
  #   }
  # }

  #Sys.setenv(ENTREZ_KEY = entrez_key)

  required_pkgs <- c("tidyverse", "readxl")
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


  library(tidyverse)
  library(readxl)
  #library(taxize)
  
  # ----- local path -----
  local <- file.path("2019_vieirasilva_naturemicrobiology_pscibd")

  # ----- File paths -----
  metadata_zip <- file.path(local,"41564_2019_483_MOESM3_ESM.xlsx.zip")
  scale_zip    <- file.path(local, "scaleandmetadata.xlsx.zip")
  qmp_zip      <- file.path(local, "QMP.matrix.tsv.zip")
  repro_counts_rds_zip <- file.path(local, "EGAS00001000001_counts.rds.zip") #REPLACE WITH THE CORRECT FILE NAME
  repro_tax_rds_zip <- file.path(local, "EGAS00001000001_taxa.rds.zip") #REPLACE WITH THE CORRECT FILE NAME

  get_rank <- function(taxon) {
    result <- classification(taxon, db = "ncbi", rows = 1)
    if (length(result) > 0 && is.data.frame(result[[1]])) {
      tax_table <- result[[1]]
      return(tax_table$rank[nrow(tax_table)])
    } else {
      return(NA)
    }
  }

  prefix_taxon <- function(original_name, rank, is_unc) {
    base_name <- sub("^unclassified_", "", original_name)
    if (!is_unc) {
      return(paste0("g_", original_name))
    } else if (is.na(rank)) {
      return(paste0("unassigned", base_name))
    } else {
      return(paste0("uc_", substr(rank, 1, 1), "_", base_name))
    }
  }
  parse_taxonomy_from_prefix <- function(prefixed_name) {
    matched <- str_match(prefixed_name, "^(g|uc_f|uc_o|uc_c|uc_p|uc_k|unassigned)_?(.*)")
    prefix <- matched[, 2]
    taxon_name <- matched[, 3]
    tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
    result <- setNames(rep(NA, length(tax_levels)), tax_levels)
    if (!is.na(prefix)) {
      rank <- tax_prefixes[[prefix]]
      if (!is.null(rank) && rank %in% names(result)) {
        result[[rank]] <- taxon_name
      }
    }
    return(result)
  }

  # ----- Initialize objects -----
  counts_original      <- NA
  proportions_original <- NA
  tax_original         <- NA
  scale       <- NA
  metadata    <- NA
  proportions_reprocessed <- NA
  counts_reprocessed <- NA
  tax_reprocessed <- NA
  proportions_reprocessed2 <- NA
  counts_reprocessed2 <- NA
  tax_reprocessed2 <- NA

  # ----- Read SCALE -----
  if (!is.na(scale_zip) && file.exists(scale_zip)) {
    temp_dir <- tempdir()
    scale_files <- unzip(scale_zip, list = TRUE)
    if (nrow(scale_files) > 0) {
      scale_xlsx <- scale_files$Name[1]
      unzip(scale_zip, files = scale_xlsx, exdir = temp_dir, overwrite = TRUE)
      scale_path <- file.path(temp_dir, scale_xlsx)
      scale_df <- read_excel(scale_path, sheet = 1)
      rename_dict <- c("Anonymised ID" = "ID")
      scale_df <- scale_df %>%
        dplyr::rename(!!!setNames(names(rename_dict), rename_dict))
      scale <- scale_df %>%
        dplyr::select(ID, `Average faecal cell count (cells/g)`) %>% 
        mutate(log2_FC_cell_g = ifelse(`Average faecal cell count (cells/g)` > 0, log2(`Average faecal cell count (cells/g)`), NA)) %>%
        mutate(log10_FC_cell_g = ifelse(`Average faecal cell count (cells/g)` > 0, log10(`Average faecal cell count (cells/g)`), NA))
      cleanup_tempfiles(temp_dir)
    }
  }

  # ----- Read METADATA -----
  if (!is.na(metadata_zip) && file.exists(metadata_zip)) {
    temp_dir <- tempdir()
    metadata_files <- unzip(metadata_zip, list = TRUE)
    metadata_xlsx <- metadata_files$Name[1]
    unzip(metadata_zip, files = metadata_xlsx, exdir = temp_dir, overwrite = TRUE)
    metadata_path <- file.path(temp_dir, metadata_xlsx)
    metadata <- read_excel(metadata_path, sheet = 2) %>%
      as.data.frame() %>%
      dplyr::rename(ID = `Anonymised ID`)
    cleanup_tempfiles(temp_dir)
  }
  

  # IF WE EVER GET THE DATA FROM THE AUTHORS ...... THEN WE CAN DELETE:collapse
  proportions_original <- read_zipped_table(file.path(local, "VieiraSilva_2019_16S.csv.zip"))
  tax_original <- data.frame(Taxa = colnames(proportions_original), stringsAsFactors = FALSE)
  if (!raw) {
    aligned = rename_and_align(proportions_original = proportions_original, metadata=metadata, scale=scale, by_col = "ID", align = align, study_name = basename(local))
    proportions_original <- aligned$proportions_original
  }
  #temp_dir <- tempdir("repro)
  #dir.create(temp_dir)

  # # ----- Read QMP and compute proportions ----- # ORIGINAL PROCESSING OF QMP MATRIX
  # if (!is.na(qmp_zip) && file.exists(qmp_zip)) {
  #   qmp_files <- unzip(qmp_zip, list = TRUE)
  #   if (nrow(qmp_files) > 0) {
  #     qmp_filename <- qmp_files$Name[1]
  #     unzip(qmp_zip, files = qmp_filename, exdir = temp_dir, overwrite = TRUE)
  #     qmp_path <- file.path(temp_dir, qmp_filename)

  #     qmp_matrix <- read.table(
  #       qmp_path,
  #       header = TRUE,
  #       sep = "\t",
  #       row.names = 1,
  #       check.names = FALSE,
  #       stringsAsFactors = FALSE
  #     )

  #     if (!is.null(scale) && !all(is.na(scale))) {
  #       common_samples <- intersect(colnames(qmp_matrix), scale$ID)
  #       scale_values <- scale %>% filter(ID %in% common_samples)
  #       scale_vector <- setNames(scale_values[["Average faecal cell count (cells/g)"]], scale_values$ID)
  #       proportions <- sweep(qmp_matrix[, names(scale_vector), drop = FALSE], 1, scale_vector, FUN = "/")
  #     } else {
  #       proportions <- qmp_matrix / rowSums(qmp_matrix)
  #     }

  #     # ----- Taxonomic rank handling -----
  #     col_names <- colnames(proportions)
  #     original_col_names <- col_names
  #     is_unclassified <- grepl("^unclassified_", col_names)
  #     clean_names <- sub("^unclassified_", "", col_names)
  #     ranks <- sapply(clean_names, get_rank)
  #     rank_mapping <- setNames(ranks, col_names)
  #     new_col_names <- mapply(prefix_taxon, original_col_names, rank_mapping, is_unclassified)
  #     manual_corrections <- list(
  #       "unassigned_Marinimicrobia" = "uc_p_Marinimicrobia",
  #       "unclassified_Blastocatella" = "uc_c_Blastocatella",
  #       "unclassified_Elioraea" = "g_Elioraea",
  #       "unclassified_Pacearchaeota" = "uc_p_Pacearchaeota",
  #       "unclassified_Bryobacter" = "g_Bryobacter",
  #       "unclassified_Geminicoccus" = "g_Geminicoccus",
  #       "unclassified_Vallitalea" = "g_Vallitalea",
  #       "unclassified_Candidatus" = "unassigned",
  #       "unclassified_Chloroplast" = "uc_c_Oxyphotobacteria",
  #       "unclassified_Natranaerovirga" = "uc_f_Natranaerovirgaceae"
  #     )

  #     for (taxon in names(manual_corrections)) {
  #       if (taxon %in% names(new_col_names)) {
  #         new_col_names[taxon] <- manual_corrections[[taxon]]
  #       }
  #     }

  #     proportions <- proportions %>%
  #       rename(!!!setNames(names(new_col_names), new_col_names))

  #     # ----- Construct taxonomy table -----
  #     tax_prefixes <- c(
  #       g = "Genus", uc_f = "Family", uc_o = "Order", uc_c = "Class",
  #       uc_p = "Phylum", uc_k = "Kingdom", unassigned = "Unassigned"
  #     )

  #     tax_rows <- lapply(names(new_col_names), function(orig_name) {
  #       prefixed <- new_col_names[orig_name]
  #       parse_taxonomy_from_prefix(prefixed)
  #     })

  #     tax <- bind_rows(tax_rows)
  #     rownames(tax) <- names(new_col_names)
  #     tax <- tibble::rownames_to_column(tax, var = "Taxon")
  #   }
  # }
  # cleanup_tempfiles(temp_dir)

  # ----- Reprocessed counts from RDS ZIP -----
  if (file.exists(repro_counts_rds_zip) && file.exists(repro_tax_rds_zip)) {
    temp_dir <- tempdir("repro")
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
    unzipped = unzip(repro_tax_rds_zip, exdir = temp_dir, overwrite = TRUE)
    tax_file <- unzipped[grep("_taxa\\.rds$", unzipped, ignore.case = TRUE)][1]
    if (is.na(tax_file)) stop("No *_taxa.rds file found after unzip")
    tax_reprocessed <- as.data.frame(readRDS(tax_file))
    
    # ----- Convert sequences to lowest rank taxonomy found and update key -----
    tax_reprocessed = make_taxa_label(tax_reprocessed)

    # ----- Convert accessions to sample IDs / Sequences to Taxa -----
    if (!raw) {
        aligned = rename_and_align(counts_reprocessed = counts_reprocessed, metadata=metadata, scale=scale, by_col = "ID", align = align, study_name = basename(local))
        counts_reprocessed <- aligned$counts_reprocessed
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

  metadata <- metadata %>% mutate(
    across(c(Diagnosis, Gender, `DMM Enterotype`), 
           ~ factor(na_if(., ""))),
    Gender = factor(Gender, 
                   levels = c("F", "M"),
                   labels = c("Female", "Male")),
    BMI = as.numeric(BMI),
    across(where(is.numeric), as.numeric)
  )

  if (!raw) {
    counts_original = fill_na_zero_numeric(counts_original)
    counts_reprocessed = fill_na_zero_numeric(counts_reprocessed)
    proportions_original = fill_na_zero_numeric(proportions_original)
    proportions_reprocessed = fill_na_zero_numeric(proportions_reprocessed)
    proportions_reprocessed2 = fill_na_zero_numeric(proportions_reprocessed2)
    counts_reprocessed2 = fill_na_zero_numeric(counts_reprocessed2)
  }

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
}
