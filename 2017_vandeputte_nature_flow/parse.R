# parse.R

require(readxl)
require(stringr)

parse_2017_vandeputte_nature_flow <- function() {
  local <- "2017_vandeputte_nature_flow/"  
  
  # Initialize the output list
  out <- list()

  #-------------------------------------------------------------
  # Read counts
  #-------------------------------------------------------------
    counts_file <- paste0(local, "OTU_nochim.csv.zip")
    if (file.exists(counts_file)) {
    message("Reading counts from ", counts_file)
    
    # Read counts file
    counts_df <- read.table(unz(counts_file, "OTU_nochim.csv"), 
                            header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    counts <- as.matrix(counts_df)
    out$counts <- counts
    
    # Calculate row sums and proportions
    row_sums <- rowSums(counts)
    proportions <- sweep(counts, 1, row_sums, FUN = "/")
    proportions[is.nan(proportions)] <- 0  # Replace NaN values with 0
    out$proportions <- proportions
    } else {
    warning("Counts file not found: ", counts_file)
    }

  #-------------------------------------------------------------
  # Read tax table
  #-------------------------------------------------------------

    taxa_file <- paste0(local, "taxa_assignments.tsv.zip")
    if (file.exists(taxa_file)) {
    message("Reading taxa assignments from ", taxa_file)
    
    # Read taxa file
    taxa_df <- read.table(unz(taxa_file, "taxa_assignments.tsv"), 
                            header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    tax <- as.matrix(taxa_df)
    out$tax <- tax
    } else {
    warning("Taxa file not found: ", taxa_file)
    }

  #-------------------------------------------------------------
  # Read scale information
  #-------------------------------------------------------------
    scale_file <- paste0(local, "cellcountstotal.zip")
    if (file.exists(scale_file)) {
        message("Reading scale from ", scale_file)
        scaledata_df <- read.table(unz(scale_file, "cellcountstotal.csv"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
        scaledata <- as.data.frame(scaledata_df, stringsAsFactors = FALSE)
        out$scale <- scaledata
    }

  #-------------------------------------------------------------
  # Read metadata
  #-------------------------------------------------------------
  
  # waiting on authors.

  #-------------------------------------------------------------
  # Checks
  #-------------------------------------------------------------
    if ("counts" %in% names(out) && "tax" %in% names(out)) {
    if (all(colnames(out$counts) %in% rownames(out$tax)) &&
        all(rownames(out$tax) %in% colnames(out$counts))) {
        message("Counts and taxa are properly aligned.")
    } else {
        warning("Mismatch detected between counts columns and taxa rows.")
    }
    }

  return(out)
}
