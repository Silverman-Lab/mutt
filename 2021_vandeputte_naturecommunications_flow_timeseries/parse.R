# parse_mydata.R

require(readxl)
require(stringr)

parse_2021_vandeputte_naturecommunications_flow_timeseries <- function() {
  local <- "2021__naturecommunications_flow_timeseries/"  
  
  # Initialize the output list
  out <- list()

  #-------------------------------------------------------------
  # Read counts
  #-------------------------------------------------------------
  counts_file <- paste0(local, "Vandeputte_2021_16S.tsv.zip")
  if (file.exists(counts_file)) {
    message("Reading counts from ", counts_file)
    counts_df <- read.table(counts_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    counts <- as.matrix(counts_df)
    out$counts <- counts
    row_sums <- rowSums(counts) 
    proportions <- sweep(counts, 1, row_sums, FUN = "/") 
    proportions[is.nan(proportions)] <- 0
    out$proportions <- proportions
  }

  #-------------------------------------------------------------
  # Read scale information
  #-------------------------------------------------------------
  scale_file <- paste0(local, "Vandeputte_2021_load.tsv.zip")
  if (file.exists(scale_file)) {
    message("Reading scale from ", scale_file)
    scaledata_df <- read.table(scale_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    scaledata <- as.data.frame(scaledata_df, stringsAsFactors = FALSE)
    out$scale <- scaledata
  }

  #-------------------------------------------------------------
  # Read tax table
  #-------------------------------------------------------------
  zipfile <- paste0(local, "41467_2021_27098_MOESM3_ESMx.xlsx.zip")
  xlsx_in_zip <- "41467_2021_27098_MOESM3_ESMx.xlsx" 
  if (!file.exists(zipfile)) {
    stop("Zip file does not exist: ", zipfile)
  }  
  tmp_dir <- tempdir()
  utils::unzip(zipfile, files = xlsx_in_zip, exdir = tmp_dir, overwrite = TRUE)
  extracted_xlsx_path <- file.path(tmp_dir, xlsx_in_zip)
  if (!file.exists(extracted_xlsx_path)) {
    stop("Failed to extract '", xlsx_in_zip, "' from '", zipfile, "'.")
  }
  library(readxl)
  sheet_name <- "S1-7"
  tax <- readxl::read_xlsx(extracted_xlsx_path, sheet = sheet_name)
  out$tax <- tax

  #-------------------------------------------------------------
  # Read metadata
  #-------------------------------------------------------------
  
  # waiting on authors.

  return(out)
}
