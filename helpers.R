# helpers.R

#' Rename and align count matrices based on metadata and scale
#'
#' This function renames the rownames of `counts_reprocessed` using the `by_col`
#' column(s) in `metadata`, and aligns both `counts_reprocessed` and `counts_original`
#' to the samples listed in `scale[[by_col]]` if `align = TRUE`.
#'
#' @param counts_reprocessed A matrix or dataframe with rownames = Accessions (optional)
#' @param counts_original A matrix or dataframe with rownames = sample identifiers (optional)
#' @param proportions_original A matrix or dataframe with rownames = sample identifiers (optional)
#' @param metadata A dataframe with columns `Accession` and the desired renaming column(s)
#' @param scale A dataframe containing the list of valid samples
#' @param by_col A character vector specifying the column(s) in both `metadata` and `scale` to align by
#' @param align Logical; if FALSE, no renaming or alignment is performed
#' @param study_name Optional name of the study for logging purposes
#'
#' @return A list with three elements: `$reprocessed`, `$original`, and `$proportions_original`
#' @export
rename_and_align <- function(counts_reprocessed = NULL,
                             counts_original = NULL,
                             proportions_original = NULL,
                             metadata,
                             scale,
                             by_col,
                             align = TRUE,
                             study_name = NULL) {
  # Validate by_col is a character vector
  stopifnot(is.character(by_col))
  stopifnot(all(by_col %in% colnames(scale)))
  
  # Create target samples based on all by_col columns
  target_samples <- if (align) {
    # Get unique combinations of all by_col values
    unique(do.call(paste, c(scale[by_col], sep = "_")))
  } else NULL

  if (!is.null(counts_reprocessed) &&
      (is.matrix(counts_reprocessed) || is.data.frame(counts_reprocessed))) {

    stopifnot("Accession" %in% colnames(metadata), all(by_col %in% colnames(metadata)))
    
    # Create mapping using all by_col columns
    acc_to_name <- setNames(
      do.call(paste, c(metadata[by_col], sep = "_")),
      metadata$Accession
    )
    
    mapped_names <- acc_to_name[rownames(counts_reprocessed)]
    valid <- !is.na(mapped_names)
    dropped <- sum(!valid)
    if (dropped > 0)
      message("\n", study_name, " | counts_reprocessed: dropped ", dropped, " unmatched")

    counts_reprocessed <- counts_reprocessed[valid, , drop = FALSE]
    mapped_names <- make.unique(mapped_names[valid])
    rownames(counts_reprocessed) <- mapped_names
    
    # Update metadata and scale with unique combinations
    metadata$combined_col <- make.unique(acc_to_name[metadata$Accession])
    scale$combined_col <- make.unique(do.call(paste, c(scale[by_col], sep = "_")))

    if (align) {
      before_n <- nrow(counts_reprocessed)
      counts_reprocessed <- counts_reprocessed[rownames(counts_reprocessed) %in% target_samples, , drop = FALSE]
      after_n <- nrow(counts_reprocessed)
      if (before_n != after_n)
        message("\n", study_name, " | counts_reprocessed: dropped ", before_n - after_n, " unaligned")
    }
  }

  if (!is.null(counts_original) &&
      (is.matrix(counts_original) || is.data.frame(counts_original)) &&
      align) {
    before_n <- nrow(counts_original)
    counts_original <- counts_original[rownames(counts_original) %in% target_samples, , drop = FALSE]
    after_n <- nrow(counts_original)
    if (before_n != after_n)
      message("\n", study_name, " | counts_original: dropped ", before_n - after_n, " unaligned")
  }

  if (!is.null(proportions_original) &&
      (is.matrix(proportions_original) || is.data.frame(proportions_original)) &&
      align) {
    before_n <- nrow(proportions_original)
    proportions_original <- proportions_original[rownames(proportions_original) %in% target_samples, , drop = FALSE]
    after_n <- nrow(proportions_original)
    if (before_n != after_n)
      message("\n", study_name, " | proportions_original: dropped ", before_n - after_n, " unaligned")
  }

  list(
    reprocessed = counts_reprocessed,
    counts_original = counts_original,
    proportions_original = proportions_original
  )
}


#' Read a single text file from a ZIP archive
#'
#' Unzips the first file inside a `.zip` archive and reads it as a data frame.
#' Defaults to treating the first column as rownames.
#'
#' @param zip_path Path to the .zip file
#' @param sep Field separator for the text file (default is comma)
#' @param header Logical; whether the file has a header row
#' @param row.names Index or name of column to use as row names
#' @param check.names Logical; whether to check and fix column names
#'
#' @return A data frame read from the zipped file, or NA if the file does not exist
#' @export
read_zipped_table <- function(zip_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE) {
  if (file.exists(zip_path)) {
    inner_file <- unzip(zip_path, list = TRUE)$Name[1]
    con <- unz(zip_path, inner_file)
    read.table(con, sep = sep, header = header, row.names = row.names,
               check.names = check.names, stringsAsFactors = FALSE)
  } else {
    warning(paste("File not found:", zip_path))
    return(NA)
  }
}

#' Generate a unique taxonomic label from taxonomic ranks
#'
#' Creates a `Taxa` column in the dataframe based on taxonomic hierarchy.
#' If `Species` column is present and not missing, it is used as the most specific label.
#' Otherwise, the function uses `Genus`, or the highest classified rank available.
#'
#' @param df A dataframe containing at least the columns: Kingdom, Phylum, Class, Order, Family, Genus.
#'           If `Species` is present, it will be used as the terminal label.
#'
#' @return The input dataframe with an added `Taxa` column
#' @export
make_taxa_label <- function(df) {
  tax_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  prefixes  <- c("k", "p", "c", "o", "f", "g")

  if (!all(tax_ranks %in% colnames(df))) {
    stop("Dataframe must contain columns: ", paste(tax_ranks, collapse = ", "))
  }

  if ("Species" %in% colnames(df)) {
    tax_ranks <- c(tax_ranks, "Species")
    prefixes  <- c(prefixes, "s")
  }

  df[tax_ranks] <- lapply(df[tax_ranks], function(x) {
    x[is.na(x) | trimws(x) == ""] <- "unclassified"
    x
  })

  df$Taxa <- apply(df[, tax_ranks], 1, function(tax_row) {
    if ("Species" %in% names(tax_row) && tax_row["Species"] != "unclassified") {
      return(paste0("s_", tax_row["Species"]))
    }
    if (tax_row["Genus"] != "unclassified") {
      return(paste0("g_", tax_row["Genus"]))
    }
    for (i in (length(tax_ranks) - 2):1) {  # skip Genus and Species
      if (tax_row[i] != "unclassified") {
        return(paste0("uc_", prefixes[i], "_", tax_row[i]))
      }
    }
    return("unclassified")
  })

  return(df)
}

#' Replace NA with 0 in numeric elements of data structures
#'
#' Recursively fills `NA` values with 0 in numeric columns of data frames, numeric matrices,
#' or nested lists containing such objects. If `x` is NULL, NA, or missing, it returns `x` unchanged.
#'
#' @param x A data frame, numeric matrix, or nested list of such objects (can be NULL or NA)
#'
#' @return The same object as `x`, with NA values in numeric components replaced by 0
#' @export
fill_na_zero_numeric <- function(x) {
  if (missing(x)) return(NULL)
  if (is.data.frame(x)) {
      x[] <- lapply(x, function(y) if (is.numeric(y)) replace(y, is.na(y), 0) else y)
  } else if (is.matrix(x) && is.numeric(x)) {
      x[is.na(x)] <- 0
  } else if (is.list(x)) {
      x <- lapply(x, fill_na_zero_numeric)
  }
  x
}

#' Read a CSV file from a ZIP archive
#'
#' Extracts and reads the first CSV file found in a ZIP archive. The function
#' assumes the first column contains row names and reads the file with those
#' as row names.
#'
#' @param zip_path Character string specifying the path to the ZIP file
#'
#' @return A data frame with row names from the first column, or NA if the file
#'         doesn't exist
#'
#' @examples
#' \dontrun{
#' # Read a CSV from a ZIP file
#' data <- read_zipped_csv("path/to/archive.zip")
#' }
#'
#' @export
read_zipped_csv <- function(zip_path) {
    if (file.exists(zip_path)) {
        csv_file <- unzip(zip_path, list = TRUE)$Name[1]
        read.csv(unz(zip_path, csv_file), row.names = 1, check.names = FALSE)
    } else {
        warning(paste("File not found:", zip_path))
        return(NA)
    }
}

#' Read an Excel file from a ZIP archive
#'
#' Extracts and reads an Excel (.xlsx) file from a ZIP archive. The function
#' can specify which sheet to read and how many rows to skip.
#'
#' @param zipfile Character string specifying the path to the ZIP file
#' @param sheet Name or index of the sheet to read (default: NULL, reads first sheet)
#' @param skip Number of rows to skip at the start of the sheet (default: 0)
#' @param tmp Directory to use for temporary file extraction (default: tempdir())
#'
#' @return A data frame containing the contents of the specified Excel sheet
#'
#' @examples
#' \dontrun{
#' # Read the first sheet of an Excel file from a ZIP archive
#' data <- read_xlsx_zip("path/to/archive.zip")
#' 
#' # Read a specific sheet, skipping the first 2 rows
#' data <- read_xlsx_zip("path/to/archive.zip", sheet = "Sheet2", skip = 2)
#' }
#'
#' @export
read_xlsx_zip <- function(zipfile,
                         sheet = NULL,
                         skip  = 0,
                         tmp   = tempdir()) {
    if (!file.exists(zipfile)) {
        stop("Zip file not found: ", zipfile)
    }
    zinfo <- utils::unzip(zipfile, list = TRUE)
    xlsx_name <- zinfo$Name[grepl("\\.xlsx$", zinfo$Name)][1]
    if (is.na(xlsx_name)) {
        stop("No .xlsx file found inside ", zipfile)
    }
    extracted <- utils::unzip(zipfile,
                            files     = xlsx_name,
                            exdir     = tmp,
                            overwrite = TRUE,
                            junkpaths = TRUE)
    xlsx_path <- file.path(tmp, basename(xlsx_name))
    if (!file.exists(xlsx_path)) {
        stop("Extraction failed—'", xlsx_path, "' does not exist")
    }
    dat <- readxl::read_xlsx(xlsx_path,
                            sheet = sheet,
                            skip  = skip)
    return(dat)
}

#' Clean up temporary files
#'
#' Safely removes temporary files and directories from the filesystem.
#' This function is useful for cleaning up after operations that create
#' temporary files, such as unzipping archives or extracting data.
#'
#' @param temp_paths Character vector of file or directory paths to remove
#'
#' @return NULL (invisibly)
#'
#' @examples
#' \dontrun{
#' # Clean up temporary files after processing
#' temp_files <- c("temp1.txt", "temp2.csv")
#' cleanup_tempfiles(temp_files)
#' }
#'
#' @export
cleanup_tempfiles <- function(temp_paths) {
    for (p in temp_paths) {
        if (file.exists(p)) {
            unlink(p, recursive = TRUE, force = TRUE)
        }
    }
}

#' Convert first row to column names
#'
#' Takes a data frame and uses its first row as column names, then removes
#' that row from the data. This is useful when reading data where the header
#' is in the first row of the data rather than as proper column names.
#'
#' @param df A data frame whose first row contains the desired column names
#'
#' @return A data frame with the first row converted to column names and removed
#'         from the data
#'
#' @examples
#' \dontrun{
#' # Convert first row to column names
#' df <- data.frame(
#'   c("Name", "John", "Jane"),
#'   c("Age", "30", "25")
#' )
#' df <- first_row_to_colnames(df)
#' # Result: data frame with columns "Name" and "Age"
#' }
#'
#' @export
first_row_to_colnames <- function(df) {
    colnames(df) <- as.character(df[1, ])
    df <- df[-1, ]
    return(df)
}

#' Build a standardized taxonomy table from various formats
#'
#' Converts taxonomy information from different classification methods (QIIME, SINTAX, UTAX, VSEARCH)
#' into a standardized format with consistent column names and formatting. The function handles
#' different input formats and cleaning steps specific to each method.
#'
#' @param df A data frame containing taxonomy information
#' @param method Character string specifying the classification method. One of:
#'   \itemize{
#'     \item "qiime" - QIIME2 sklearn classifier output
#'     \item "sintax" - SINTAX output with 80% confidence
#'     \item "sintax_full" - Full SINTAX output
#'     \item "utax" - UTAX classifier output
#'     \item "vsearch" - VSEARCH usearch_global output
#'   }
#'
#' @return A data frame with standardized taxonomy columns:
#'   \itemize{
#'     \item OTU_ID - Original sequence identifier
#'     \item Kingdom - Kingdom classification
#'     \item Phylum - Phylum classification
#'     \item Class - Class classification
#'     \item Order - Order classification
#'     \item Family - Family classification
#'     \item Genus - Genus classification
#'     \item Species - Species classification
#'   }
#'
#' @examples
#' \dontrun{
#' # Process QIIME2 taxonomy
#' tax_table <- build_taxonomy_table(my_data, method = "qiime")
#' 
#' # Process SINTAX taxonomy
#' tax_table <- build_taxonomy_table(my_data, method = "sintax")
#' }
#'
#' @export
build_taxonomy_table <- function(df, method = c("qiime", "sintax", "sintax_full", "utax", "vsearch")) {
    method <- match.arg(method)
    if (method == "qiime") {
      tax_source <- df %>%
        select(OTU_ID, taxonomy = qiime_sklearn)
      if (any(grepl("^D_\\d+__", tax_source$taxonomy, perl = TRUE))) {
        tax_df <- tax_source %>%
          separate(taxonomy,
                  into = c("D0", "D1", "D2", "D3", "D4", "D5", "D6"),
                  sep = ";", fill = "right") %>%
          transmute(
            OTU_ID,
            Kingdom = gsub("^D_0__", "", D0),
            Phylum  = gsub("^D_1__", "", D1),
            Class   = gsub("^D_2__", "", D2),
            Order   = gsub("^D_3__", "", D3),
            Family  = gsub("^D_4__", "", D4),
            Genus   = gsub("^D_5__", "", D5),
            Species = gsub("^D_6__", "", D6)
          )
      } else {
        tax_df <- tax_source %>%
          separate(taxonomy,
                  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                  sep = ";", fill = "right") %>%
          mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))
      }
      tax_df <- tax_df %>%
        mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), NA, .)))
    } else if (method == "sintax") {
      tax_df <- df %>%
        select(OTU_ID, taxonomy = `usearch_sintax_80%`) %>%
        separate(taxonomy,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ",", fill = "right") %>%
        mutate(across(Kingdom:Species, ~ gsub("^[a-z]:", "", .)))
    } else if (method == "sintax_full") {
      tax_df <- df %>%
        select(OTU_ID, taxonomy = usearch_sintax) %>%
        mutate(taxonomy = gsub("\\([^)]+\\)", "", taxonomy)) %>%
        separate(taxonomy,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ",", fill = "right") %>%
        mutate(across(Kingdom:Species, ~ gsub("^[a-z]:", "", .)))
    } else if (method == "utax") {
      tax_df <- df %>%
        select(OTU_ID, taxonomy = usearch_utax) %>%
        mutate(taxonomy = gsub(".*\\|refs\\|", "", taxonomy)) %>%
        separate(taxonomy,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";", fill = "right") %>%
        mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))
    } else if (method == "vsearch") {
      tax_df <- df %>%
        select(OTU_ID, taxonomy = vsearch_usearchglobal) %>%
        mutate(taxonomy = gsub(".*\\|refs\\|", "", taxonomy)) %>%
        separate(taxonomy,
                into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                sep = ";", fill = "right") %>%
        mutate(across(Kingdom:Species, ~ gsub("^[a-z]__", "", .)))
    }
    tax_df <- tax_df %>%
      mutate(across(Kingdom:Species, ~ ifelse(. == "" | is.na(.), NA, .)))
    return(tax_df)
}

#' Validate the structure of a dataset output
#'
#' Checks if a dataset output contains all required elements and returns a structured
#' object showing what elements are present and their types.
#'
#' @param output A list containing the dataset output
#' @param study_name Optional name of the study for error messages
#'
#' @return A list with the same structure as the input, but with element types indicated
#'
#' @examples
#' \dontrun{
#' # Validate a dataset output
#' result <- validate_output_structure(my_output, study_name = "MyStudy")
#' }
#'
#' @export
validate_output_structure <- function(output, study_name = NULL) {
    # Check if output is a list
    if (!is.list(output)) {
        stop(sprintf("%s: Output must be a list", ifelse(is.null(study_name), "", study_name)))
    }

    # Required top-level elements
    required_elements <- c("counts", "proportions", "tax", "scale", "metadata")
    missing_elements <- required_elements[!required_elements %in% names(output)]
    
    # Print warning for missing elements if study_name is provided
    if (length(missing_elements) > 0 && !is.null(study_name)) {
        warning(sprintf("%s: Missing required elements: %s", 
                       study_name,
                       paste(missing_elements, collapse = ", ")))
    }

    # Create a structured output showing what's present
    result <- list()
    for (elem in names(output)) {
        if (is.list(output[[elem]])) {
            if (all(c("original", "reprocessed") %in% names(output[[elem]]))) {
                result[[elem]] <- "LIST={original,reprocessed}"
            } else {
                result[[elem]] <- "LIST"
            }
        } else if (is.data.frame(output[[elem]])) {
            result[[elem]] <- "DATAFRAME"
        } else if (is.matrix(output[[elem]])) {
            result[[elem]] <- "MATRIX"
        } else {
            result[[elem]] <- class(output[[elem]])[1]
        }
    }

    return(result)
}

