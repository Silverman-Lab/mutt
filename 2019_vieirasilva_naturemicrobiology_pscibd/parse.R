parse_2019_vieirasilva_naturemicrobiology_pscibd <- function(raw=FALSE,entrez_key = NULL) {
  if (is.null(entrez_key) || !nzchar(entrez_key)) {
    stop(
      "You must provide a valid NCBI Entrez API key via the `entrez_key` argument.\n",
      "To create or retrieve your key:\n",
      "  1) Visit: https://www.ncbi.nlm.nih.gov/account/\n",
      "  2) Sign in or create an account\n",
      "  3) Go to 'Settings' → 'API Keys' and generate a new key\n",
      "  4) Provide that key to this function, e.g.:\n",
      "     parse_2019_viera-silva_naturemicrobiology_pscibd(entrez_key = \"YOUR_KEY_HERE\")"
    )
  }

  Sys.setenv(ENTREZ_KEY = entrez_key)

  required_pkgs <- c("tidyverse", "readxl", "taxize")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }

  library(tidyverse)
  library(readxl)
  library(taxize)
  
  # ----- local path -----
  
  local <- file.path("2019_vieirasilva_naturemicrobiology_pscibd")

  # ----- File paths -----
  metadata_zip <- file.path(local,"41564_2019_483_MOESM3_ESM.xlsx.zip")
  scale_zip    <- file.path(local, "scaleandmetadata.xlsx.zip")
  qmp_zip      <- file.path(local, "QMP.matrix.tsv.zip")

  # ----- Initialize objects -----
  counts      <- NA
  proportions <- NA
  tax         <- NA
  scale       <- NA
  metadata    <- NA

  temp_dir <- tempdir()

  # ----- Read SCALE -----
  if (!is.na(scale_zip) && file.exists(scale_zip)) {
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
        dplyr::select(ID, `Average faecal cell count (cells/g)`)
    }
  }

  # ----- Read QMP and compute proportions -----
  if (!is.na(qmp_zip) && file.exists(qmp_zip)) {
    qmp_files <- unzip(qmp_zip, list = TRUE)
    if (nrow(qmp_files) > 0) {
      qmp_filename <- qmp_files$Name[1]
      unzip(qmp_zip, files = qmp_filename, exdir = temp_dir, overwrite = TRUE)
      qmp_path <- file.path(temp_dir, qmp_filename)

      qmp_matrix <- read.table(
        qmp_path,
        header = TRUE,
        sep = "\t",
        row.names = 1,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )

      if (!is.null(scale) && !all(is.na(scale))) {
        common_samples <- intersect(colnames(qmp_matrix), scale$ID)
        scale_values <- scale %>% filter(ID %in% common_samples)
        scale_vector <- setNames(scale_values[["Average faecal cell count (cells/g)"]], scale_values$ID)
        proportions <- sweep(qmp_matrix[, names(scale_vector), drop = FALSE], 2, scale_vector, FUN = "/")
      } else {
        proportions <- qmp_matrix / rowSums(qmp_matrix)
      }

      # ----- Taxonomic rank handling -----
      col_names <- colnames(proportions)
      original_col_names <- col_names
      is_unclassified <- grepl("^unclassified_", col_names)
      clean_names <- sub("^unclassified_", "", col_names)

      get_rank <- function(taxon) {
        result <- classification(taxon, db = "ncbi", rows = 1)
        if (length(result) > 0 && is.data.frame(result[[1]])) {
          tax_table <- result[[1]]
          return(tax_table$rank[nrow(tax_table)])
        } else {
          return(NA)
        }
      }

      ranks <- sapply(clean_names, get_rank)
      rank_mapping <- setNames(ranks, col_names)

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

      new_col_names <- mapply(prefix_taxon, original_col_names, rank_mapping, is_unclassified)

      manual_corrections <- list(
        "unassigned_Marinimicrobia" = "uc_p_Marinimicrobia",
        "unclassified_Blastocatella" = "uc_c_Blastocatella",
        "unclassified_Elioraea" = "g_Elioraea",
        "unclassified_Pacearchaeota" = "uc_p_Pacearchaeota",
        "unclassified_Bryobacter" = "g_Bryobacter",
        "unclassified_Geminicoccus" = "g_Geminicoccus",
        "unclassified_Vallitalea" = "g_Vallitalea",
        "unclassified_Candidatus" = "unassigned",
        "unclassified_Chloroplast" = "uc_c_Oxyphotobacteria",
        "unclassified_Natranaerovirga" = "uc_f_Natranaerovirgaceae"
      )

      for (taxon in names(manual_corrections)) {
        if (taxon %in% names(new_col_names)) {
          new_col_names[taxon] <- manual_corrections[[taxon]]
        }
      }

      proportions <- proportions %>%
        rename(!!!setNames(names(new_col_names), new_col_names))

      # ----- Construct taxonomy table -----
      tax_prefixes <- c(
        g = "Genus", uc_f = "Family", uc_o = "Order", uc_c = "Class",
        uc_p = "Phylum", uc_k = "Kingdom", unassigned = "Unassigned"
      )

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

      tax_rows <- lapply(names(new_col_names), function(orig_name) {
        prefixed <- new_col_names[orig_name]
        parse_taxonomy_from_prefix(prefixed)
      })

      tax <- bind_rows(tax_rows)
      rownames(tax) <- names(new_col_names)
      tax <- tibble::rownames_to_column(tax, var = "Taxon")
    }
  }

  # ----- Read METADATA -----
  if (!is.na(metadata_zip) && file.exists(metadata_zip)) {
    metadata_csv <- unzip(metadata_zip, list = TRUE)$Name[1]
    metadata_con <- unz(metadata_zip, metadata_csv)
    metadata <- read.csv(metadata_con, row.names = "Anonymised ID") %>%
      as.data.frame()
  }

  #   # ----- Reprocessed counts from RDS ZIP -----
  # temp_rds <- tempfile(fileext = ".rds")
  # unzip(repro_counts_rds_zip, exdir = dirname(temp_rds), overwrite = TRUE)

  # rds_files <- list.files(dirname(temp_rds), pattern = "_counts\\.rds$", full.names = TRUE)
  # if (length(rds_files) == 0) stop("No *_counts.rds file found after unzip")
  # counts_reprocessed <- as.data.frame(readRDS(rds_files[1]))

  # # ----- Taxonomy reprocessed -----
  # temp_tax <- tempfile(fileext = ".rds")
  # unzip(repro_tax_zip, exdir = dirname(temp_tax), overwrite = TRUE)

  # tax_files <- list.files(dirname(temp_tax), pattern = "_taxa\\.rds$", full.names = TRUE)
  # if (length(tax_files) == 0) stop("No *_taxa.rds file found after unzip")
  # tax_reprocessed <- as.data.frame(readRDS(tax_files[1]))

  
  # # ----- Convert sequences to lowest rank taxonomy found and update key -----
  # make_taxa_label <- function(df) {
  #     tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  #     prefixes  <- c("k", "p", "c", "o", "f", "g")
  #     if (!all(tax_ranks %in% colnames(df))) {
  #         stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
  #     }
  #     df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
  #         x[is.na(x) | trimws(x) == ""] <- "unclassified"
  #         x
  #     })
  #     df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
  #         if (tax_row["Genus"] != "unclassified") {
  #         return(paste0("g_", tax_row["Genus"]))
  #         }
  #         for (i in (length(tax_ranks)-1):1) {  # skip Genus
  #         if (tax_row[i] != "unclassified") {
  #             return(paste0("uc_", prefixes[i], "_", tax_row[i]))
  #         }
  #         }
  #         return("unclassified")
  #     })
  #     return(df)
  # }
  # tax_reprocessed = make_taxa_label(tax_reprocessed)

  # # ----- Convert accessions to sample IDs / Sequences to Taxa -----
  # # accessions to sampleIDs is study specific: IF NEED BE

  # # taxa
  # if (!raw) {
  #     matched_taxa <- tax_reprocessed$Taxa[match(colnames(counts_reprocessed), rownames(tax_reprocessed))]
  #     colnames(counts_reprocessed) <- matched_taxa
  #     counts_reprocessed <- as.data.frame(t(rowsum(t(counts_reprocessed), group = colnames(counts_reprocessed))))
  # }

  # # proportions reprocessed
  # proportions_reprocessed = counts_reprocessed
  # proportions_reprocessed[-1] <- lapply(
  #     counts_reprocessed[-1],
  #     function(col) col / sum(col)
  # )

  return(list(
    counts = list(
      original = NA,
      reprocessed = NA
    ),
    proportions = list(
      original = proportions,
      reprocessed = NA
    ),
    tax = list(
      original = tax,
      reprocessed = NA
    ),
    scale = scale,
    metadata = metadata
  ))
}
