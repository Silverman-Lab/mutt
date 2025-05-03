parse_2024_nishijima_cell_galaxy <- function(bind_all = FALSE) {
  required_pkgs <- c("stringr", "tidyverse")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "),
         ". Please install them before running this function.")
  }
<<<<<<< Updated upstream
<<<<<<< Updated upstream

  library(stringr)
  library(tidyverse)

  # ----- Local base directory -----
  local <- file.path("2024_nishijima_cell_galaxy")

  # ----- Dataset prefixes -----
  prefixes <- c("GALA-ALD", "GALA-HP", "AlcoChallenge", "GALA-RIF", "GastricBypass",
                "HCO", "GALA-POSTBIO", "HOLABEK")

  # ----- File IDs corresponding to each prefix -----
  project_ids <- c(
    "PRJEB76661", "PRJEB76664", "PRJEB76662", "PRJEB76667",
    "PRJEB76666", "PRJEB81698", "PRJEB76668", "PRJEB81697"
  )

  # ----- File paths -----
  motus_files <- setNames(
    file.path(local, paste0(project_ids, "_motus_merged.tsv.zip")),
    prefixes
  )

  metaphlan_files <- setNames(
    file.path(local, paste0(project_ids, "_MetaPhlAn_merged.tsv.zip")),
    prefixes
  )
=======
=======
>>>>>>> Stashed changes

  library(stringr)
  library(tidyverse)

  local <- "2024_nishijima_cell_galaxy/"
  prefixes <- c("GALA-ALD", "GALA-HP", "AlcoChallenge", "GALA-RIF", "GastricBypass",
                "HCO", "GALA-POSTBIO", "HOLABEK")

  motus_files <- setNames(paste0(local, c(
    "PRJEB76661", "PRJEB76664", "PRJEB76662", "PRJEB76667",
    "PRJEB76666", "PRJEB81698", "PRJEB76668", "PRJEB81697"
  ), "_motus_merged.tsv.zip"), prefixes)

  metaphlan_files <- setNames(paste0(local, c(
    "PRJEB76661", "PRJEB76664", "PRJEB76662", "PRJEB76667",
    "PRJEB76666", "PRJEB81698", "PRJEB76668", "PRJEB81697"
  ), "_MetaPhlAn_merged.tsv.zip"), prefixes)
<<<<<<< Updated upstream
>>>>>>> Stashed changes
=======
>>>>>>> Stashed changes

  # ----- Load Scale & Metadata -----
  dat <- read.table(paste0(local, "GALAXY_load.tsv.gz"), header = TRUE)
  scale <- data.frame(count = dat$count)
  rownames(scale) <- dat$ID
  scale <- as.matrix(scale)

  metadata <- data.frame(cohort = dat$cohort)
  rownames(metadata) <- dat$ID

  # ----- Load Original Proportions Table -----
  dat <- read.table(paste0(local, "GALAXY_mOTUs_v25.tsv.gz"), header = FALSE, skip = 1)
  proportions_raw <- dat
  rownames(proportions_raw) <- dat$V1
  proportions_raw <- proportions_raw[, -1]
  colnames(proportions_raw)[colnames(proportions_raw) == "-1"] <- "unclassified"

  line1 <- readLines(paste0(local, "GALAXY_mOTUs_v25.tsv.gz"), n = 1)
  taxnames <- strsplit(line1, "\t")[[1]]
  seqids <- str_extract(taxnames[-1], "(?<=mOTU_v25_)[0-9]+")
  proportions <- t(as.matrix(proportions_raw))
  rownames(proportions) <- seqids
  colnames(proportions) <- rownames(proportions_raw)

  # ----- Original Taxonomy from Mapping -----
  tax_map <- read_tsv(paste0(local, "motus2GTDB.txt.gz"), col_names = c("mOTU_ID", "taxonomy"))
  tax_table <- tax_map %>%
    mutate(mOTU_ID = as.character(mOTU_ID)) %>%
    filter(mOTU_ID %in% rownames(proportions)) %>%
    mutate(taxonomy = str_replace_all(taxonomy, "\\|", ";")) %>%
    separate(taxonomy,
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
             sep = ";", fill = "right", extra = "drop") %>%
    column_to_rownames("mOTU_ID")

  # ----- Helper to Load Reprocessed -----
  load_reprocessed <- function(zip_path) {
    if (!file.exists(zip_path)) return(NULL)
    file_list <- unzip(zip_path, list = TRUE)
    file_name <- file_list$Name[grepl("\\.tsv$", file_list$Name)][1]
    if (is.na(file_name)) return(NULL)

    temp_dir <- tempdir()
    unzip(zip_path, files = file_name, exdir = temp_dir, overwrite = TRUE)
    df <- read_tsv(file.path(temp_dir, file_name))
    rownames(df) <- df[[1]]
    df[[1]] <- NULL
    prop <- apply(df, 2, function(col) col / sum(col))
    tax <- data.frame(taxa = rownames(df)) %>%
      mutate(taxa = str_trim(taxa)) %>%
      separate(taxa,
               into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
               sep = "\\s*;\\s*", extra = "drop", fill = "right")
    rownames(tax) <- rownames(df)
    list(counts = df, proportions = prop, tax = tax)
  }

  # ----- Load All Reprocessed Files -----
  mOTU3_reprocessed <- list()
  MetaPhlAn4_reprocessed <- list()

  for (prefix in prefixes) {
    motus_res <- load_reprocessed(motus_files[prefix])
    if (!is.null(motus_res)) {
      mOTU3_reprocessed[[prefix]] <- motus_res
    }
    metaphlan_res <- load_reprocessed(metaphlan_files[prefix])
    if (!is.null(metaphlan_res)) {
      MetaPhlAn4_reprocessed[[prefix]] <- metaphlan_res
    }
  }

  # ----- Return Based on `bind_all` Option -----
  return(list(
    counts = list(
      original = NA,
      reprocessed = list( 
          mOTU3 = if (bind_all) bind_rows(lapply(mOTU3_reprocessed, `[[`, "counts"), .id = "prefix") else lapply(mOTU3_reprocessed, `[[`, "counts"),
          MetaPhlAn4 = if (bind_all) bind_rows(lapply(MetaPhlAn4_reprocessed, `[[`, "counts"), .id = "prefix") else lapply(MetaPhlAn4_reprocessed, `[[`, "counts")
      )
    ),
    proportions = list(
      original = proportions,
      reprocessed = list(
          mOTU3 = if (bind_all) bind_rows(lapply(mOTU3_reprocessed, `[[`, "proportions"), .id = "prefix") else lapply(mOTU3_reprocessed, `[[`, "proportions"),
          MetaPhlAn4 = if (bind_all) bind_rows(lapply(MetaPhlAn4_reprocessed, `[[`, "proportions"), .id = "prefix") else lapply(MetaPhlAn4_reprocessed, `[[`, "proportions")
      )
    ),
    tax = list(
      original = tax_table,
      reprocessed = list(
          mOTU3 = if (bind_all) bind_rows(lapply(mOTU3_reprocessed, `[[`, "tax"), .id = "prefix") else lapply(mOTU3_reprocessed, `[[`, "tax"),
          MetaPhlAn4 = if (bind_all) bind_rows(lapply(MetaPhlAn4_reprocessed, `[[`, "tax"), .id = "prefix") else lapply(MetaPhlAn4_reprocessed, `[[`, "tax")
      )
    ),
    scale = scale,
    metadata = metadata
  ))
}
