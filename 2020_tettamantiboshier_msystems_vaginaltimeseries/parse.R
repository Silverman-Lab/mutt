# parse.R

require(readxl)
require(stringr)

parse_2020_tettamantiboshier_msystems_vaginaltimeseries <- function() {
  local <- "2020_tettamantiboshier_msystems_vaginaltimeseries/"  

  out <- list()

  #-------------------------------------------------------------
  # Read counts
  #-------------------------------------------------------------
  counts_file <- paste0(local, "mSystems.00777-19-st004.zip")
    if (file.exists(counts_file)) {
    message("Reading counts from ", counts_file)
    
    # Read counts file
    counts_df <- read.table(unz(counts_file, "mSystems.00777-19-st004.csv"), 
                            header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    counts <- as.matrix(counts_df)
    out$counts <- counts
    
    # Calculate row sums and proportions
    row_sums <- rowSums(counts)
    proportions <- sweep(counts, 1, row_sums, FUN = "/")
    proportions[is.nan(proportions)] <- 0  
    out$proportions <- proportions

  #-------------------------------------------------------------
  # Taxa data
  #-------------------------------------------------------------
    out$tax <- colnames(counts_df)[3:ncol(counts_df)]

  #-------------------------------------------------------------
  # Scale data
  #-------------------------------------------------------------
    scale_file <- paste0(local, "mSystems.00777-19-st002.zip")
    if (file.exists(scale_file)) {
    message("Reading counts from ", scale_file)
    
    # Read counts file
    scale_df <- read.table(unz(scale_file, "mSystems.00777-19-st002.csv"), 
                            header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)  
    out$scale <- scale_df %>%
            dplyr::select(Participant, Hours_In_Study, Total_16S)


  #-------------------------------------------------------------
  # Read metadata 
  #-------------------------------------------------------------

  
  #-------------------------------------------------------------
  # Checks for alignment
  #-------------------------------------------------------------
  
  #-------------------------------------------------------------
  # Return the combined dataset
  #-------------------------------------------------------------
  return(out)
};