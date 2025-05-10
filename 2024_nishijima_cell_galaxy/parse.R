parse_2024_nishijima_cell_galaxy <- function(raw = FALSE, align = FALSE, bind_all = FALSE) {
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

  sra_files <- c(
    "PRJEB76661" = file.path(local,"SraRunTable (40).csv.zip"),
    "PRJEB76664" = file.path(local,"SraRunTable (41).csv.zip"),
    "PRJEB76662" = file.path(local,"SraRunTable (42).csv.zip"),
    "PRJEB76667" = file.path(local,"SraRunTable (43).csv.zip"),
    "PRJEB76666" = file.path(local,"SraRunTable (44).csv.zip"),
    "PRJEB81698" = file.path(local,"SraRunTable (45).csv.zip"),
    "PRJEB76668" = file.path(local,"SraRunTable (46).csv.zip"),
    "PRJEB81697" = file.path(local,"SraRunTable (47).csv.zip")
  )

  counts_original = NA
  mOTU3_reprocessed = NA
  MetaPhlAn4_reprocessed = NA
  proportions_original = NA
  tax_original = NA
  scale = NA
  metadata = NA

  # ----- Load Scale & Metadata -----
  dat <- read_zipped_table(file.path(local, "GALAXY_load.tsv.zip"), sep="\t", header = TRUE, row.names = NULL)
  scale <- data.frame(ID = dat$ID, count = dat$count)
  rownames(scale) <- dat$ID

  metadata <- data.frame(ID = dat$ID, cohort = dat$cohort)
  rownames(metadata) <- dat$ID

  scale <- scale %>% mutate(log2_FC = ifelse(count >0, log2(count), NA)) %>% mutate(log10_FC = ifelse(count> 0, log10(count), NA))

  if (any(file.exists(sra_files))) {
    metadata <- metadata %>%
      left_join(sra_files, by = c("ID" = "Sample Name")) %>% rename(Accession = Run)
  }

  # ----- Load Original Proportions Table -----
  # First get the header line for taxnames
  con <- unz(file.path(local, "GALAXY_mOTUs_v25.tsv.zip"), unzip(file.path(local, "GALAXY_mOTUs_v25.tsv.zip"), list = TRUE)$Name[1])
  line1 <- readLines(con, n = 1)
  close(con)
  taxnames <- strsplit(line1, "\t")[[1]]
  seqids <- str_extract(taxnames[-1], "(?<=mOTU_v25_)[0-9]+")
  proportions_original <- read_zipped_table(file.path(local, "GALAXY_mOTUs_v25.tsv.zip"), sep="\t", header = FALSE, skip = 1)
  proportions_original = proportions_original %>% rownames_to_column(var = "ID")
  if (!raw) {
    aligned = rename_and_align(proportions_original = proportions_original, metadata=metadata, scale=scale, by_col="ID", align = align, study_name=basename(local))
    proportions_original = aligned$proportions_original
  }
  rownames(proportions_original) <- proportions_original$ID
  proportions_original<- proportions_original[, -1]

  # Map seqids to column names
  colnames(proportions_original) <- taxnames[-1]

  # ----- Original Taxonomy from Mapping -----
  tax_original <- read_zipped_table(file.path(local, "motus2GTDB.txt.zip"), sep="\t")

  # # ----- Helper to Load Reprocessed -----
  # load_reprocessed <- function(zip_path) {
  #   if (!file.exists(zip_path)) return(NULL)
  #   file_list <- unzip(zip_path, list = TRUE)
  #   file_name <- file_list$Name[grepl("\\.tsv$", file_list$Name)][1]
  #   if (is.na(file_name)) return(NULL)

  #   temp_dir <- tempdir()
  #   unzip(zip_path, files = file_name, exdir = temp_dir, overwrite = TRUE)
    
  #   # Ensure first column is character for row names
  #   df[[1]] <- as.character(df[[1]])
  #   rownames(df) <- df[[1]]
  #   df[[1]] <- NULL
    
  #   # Set row names
  #   rownames(df) <- row_names
    
  #   if (!raw) {
  #     if (any(file.exists(sra_files))) {
  #       aligned = rename_and_align(counts_reprocessed = df, metadata=metadata, scale=scale, by_col="ID", align = align, study_name=basename(local))
  #       df = aligned$reprocessed
  #     } else {
  #       aligned = rename_and_align(counts_original = df, metadata=metadata, scale=scale, by_col="ID", align = align, study_name=basename(local))
  #       df = aligned$counts_original
  #     }
  #   }
    
  #   # Calculate proportions using matrix operations
  #   row_sums <- rowSums(df, na.rm = TRUE)
  #   prop <- sweep(df, 1, row_sums, FUN = "/")
  #   prop[is.nan(prop)] <- 0
    
  #   # Handle taxonomy
  #   tax <- data.frame(taxa = rownames(df), stringsAsFactors = FALSE) %>%
  #     mutate(taxa = str_trim(taxa)) %>%
  #     separate(taxa,
  #              into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"),
  #              sep = "\\s*;\\s*", extra = "drop", fill = "right")
  #   rownames(tax) <- rownames(df)
    
  #   list(counts = df, proportions = prop, tax = tax)
  # }

  # # ----- Load All Reprocessed Files -----
  # mOTU3_reprocessed <- list()
  # MetaPhlAn4_reprocessed <- list()

  # if (any(file.exists(motus_files)) || any(file.exists(metaphlan_files))) {
  #   for (prefix in prefixes) {
  #     motus_res <- load_reprocessed(motus_files[prefix])
  #     if (!is.null(motus_res)) {
  #       mOTU3_reprocessed[[prefix]] <- motus_res
  #     }
  #     metaphlan_res <- load_reprocessed(metaphlan_files[prefix])
  #     if (!is.null(metaphlan_res)) {
  #       MetaPhlAn4_reprocessed[[prefix]] <- metaphlan_res
  #     }
  #   }
  # }

  # ----- Return Based on `bind_all` Option -----
  return(list(
    counts = list(
      original = counts_original,
      reprocessed = list( 
          mOTU3 = if (!is.na(mOTU3_reprocessed)[1]) {
              if (length(mOTU3_reprocessed) > 0) {
                  if (bind_all) {
                      do.call(rbind, lapply(mOTU3_reprocessed, function(x) {
                          df <- as.data.frame(x$counts)
                          df$prefix <- names(x)[1]
                          df
                      }))
                  } else {
                      lapply(mOTU3_reprocessed, `[[`, "counts")
                  }
              } else {
                  NA
              }
          } else {
              NA
          },
          MetaPhlAn4 = if (!is.na(MetaPhlAn4_reprocessed)[1]) {
              if (length(MetaPhlAn4_reprocessed) > 0) {
                  if (bind_all) {
                      do.call(rbind, lapply(MetaPhlAn4_reprocessed, function(x) {
                          df <- as.data.frame(x$counts)
                          df$prefix <- names(x)[1]
                          df
                      }))
                  } else {
                      lapply(MetaPhlAn4_reprocessed, `[[`, "counts")
                  }
              } else {
                  NA
              }
          } else {
              NA
          }
      )
    ),
    proportions = list(
      original = proportions_original,
      reprocessed = list(
          mOTU3 = if (!is.na(mOTU3_reprocessed)[1]) {
              if (length(mOTU3_reprocessed) > 0) {
                  if (bind_all) {
                      do.call(rbind, lapply(mOTU3_reprocessed, function(x) {
                          df <- as.data.frame(x$proportions)
                          df$prefix <- names(x)[1]
                          df
                      }))
                  } else {
                      lapply(mOTU3_reprocessed, `[[`, "proportions")
                  }
              } else {
                  NA
              }
          } else {
              NA
          },
          MetaPhlAn4 = if (!is.na(MetaPhlAn4_reprocessed)[1]) {
              if (length(MetaPhlAn4_reprocessed) > 0) {
                  if (bind_all) {
                      do.call(rbind, lapply(MetaPhlAn4_reprocessed, function(x) {
                          df <- as.data.frame(x$proportions)
                          df$prefix <- names(x)[1]
                          df
                      }))
                  } else {
                      lapply(MetaPhlAn4_reprocessed, `[[`, "proportions")
                  }
              } else {
                  NA
              }
          } else {
              NA
          }
      )
    ),
    tax = list(
      original = tax_original,
      reprocessed = list(
          mOTU3 = if (!is.na(mOTU3_reprocessed)[1]) {
              if (length(mOTU3_reprocessed) > 0) {
                  if (bind_all) {
                      do.call(rbind, lapply(mOTU3_reprocessed, function(x) {
                          df <- as.data.frame(x$tax)
                          df$prefix <- names(x)[1]
                          df
                      }))
                  } else {
                      lapply(mOTU3_reprocessed, `[[`, "tax")
                  }
              } else {
                  NA
              }
          } else {
              NA
          },
          MetaPhlAn4 = if (!is.na(MetaPhlAn4_reprocessed)[1]) {
              if (length(MetaPhlAn4_reprocessed) > 0) {
                  if (bind_all) {
                      do.call(rbind, lapply(MetaPhlAn4_reprocessed, function(x) {
                          df <- as.data.frame(x$tax)
                          df$prefix <- names(x)[1]
                          df
                      }))
                  } else {
                      lapply(MetaPhlAn4_reprocessed, `[[`, "tax")
                  }
              } else {
                  NA
              }
          } else {
              NA
          }
      )
    ),
    scale = scale,
    metadata = metadata
  ))
}
