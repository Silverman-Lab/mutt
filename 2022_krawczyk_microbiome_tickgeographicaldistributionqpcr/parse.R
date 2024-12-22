# parse.R

require(readxl)
require(stringr)

parse_2022_krawczyk_microbiome_tickgeographicaldistributionqpcr <- function() {
  local <- "2022_krawczyk_microbiome_tickgeographicaldistributionqpcr/"  

  out <- list()

  #-------------------------------------------------------------
  # Read counts
  #-------------------------------------------------------------

  countsfile <- paste0(local, "40168_2022_1276_MOESM3_ESM.zip")
  if (file.exists(countsfile)) {
    message("Extracting and reading data from ", countsfile)

    temp_dir <- tempdir()
    unzip(countsfile, exdir = temp_dir)
    excel_files <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)
    
    if (length(excel_files) > 0) {
      counts_file <- excel_files[1]  
      counts_sheets <- c("TableS4", "TableS5")
      
      counts_list <- lapply(counts_sheets, function(sheet) {
        read_excel(counts_file, sheet = sheet)
      })
      names(counts_list) <- counts_sheets

        if ("TableS4" %in% names(counts_list)) {
        excluded_columns <- c("Taxonomy")  
        out$counts <- counts_list$"TableS4" %>%
            dplyr::select(-all_of(excluded_columns))
  #-------------------------------------------------------------
  # Taxa data
  #-------------------------------------------------------------
        out$tax <- counts_list$"TableS4" %>%
            dplyr::select(Name, Taxonomy)
    
        } else {
        warning("TableS4 not found in the Excel file.")
        }
  #-------------------------------------------------------------
  # Scale data
  #-------------------------------------------------------------
        # Add TableS5 directly (if needed in its entirety)
        if ("TableS5" %in% names(counts_list)) {
        out$scale <- counts_list$"TableS5"
        } else {
        warning("TableS5 not found in the Excel file.")
        }

      } else {
      warning("No .xlsx file found in ", metadata_zip)
    }
  } else {
    warning("Metadata zip file not found: ", metadata_zip)
  }


  #-------------------------------------------------------------
  # Read metadata 
  #-------------------------------------------------------------
  metadata_zip <- paste0(local, "40168_2022_1276_MOESM1_ESM.zip")
  if (file.exists(metadata_zip)) {
    message("Extracting and reading metadata from ", metadata_zip)

    temp_dir <- tempdir()
    unzip(metadata_zip, exdir = temp_dir)
    excel_files <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)
    
    if (length(excel_files) > 0) {
      metadata_file <- excel_files[1]  
      metadata_sheets <- c("Table S1", "Table S3")
      
      metadata_list <- lapply(metadata_sheets, function(sheet) {
        read_excel(metadata_file, sheet = sheet)
      })
      names(metadata_list) <- metadata_sheets
      out$metadata <- metadata_list
    } else {
      warning("No .xlsx file found in ", metadata_zip)
    }
  } else {
    warning("Metadata zip file not found: ", metadata_zip)
  }
  
  #-------------------------------------------------------------
  # Checks for alignment
  #-------------------------------------------------------------
  
  #-------------------------------------------------------------
  # Return the combined dataset
  #-------------------------------------------------------------
  return(out)
}