# parse.R

require(openxlsx)
require(stringr)

parse_2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification <- function() {
  local <- "2021_rao_nature_mkspikeseqmetagenomicmultiplescalequantification/"  
  
    # Initialize the output list
    out <- list()

    #-------------------------------------------------------------
    # Read counts
    #-------------------------------------------------------------

    #*** Need to fix for archaea ***

    counts_file <- paste0(local, "combined_p1_p2.xlsx.zip")

    if (file.exists(counts_file)) {
    message("Reading proportions from ", counts_file)
    
    temp_dir <- tempdir()
    unzip(counts_file, exdir = temp_dir)
    excel_file <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)
    
    sheets <- c("fungi", "bacteria", "archaea", "archaea_side")
    data_list <- lapply(sheets, function(sheet) {
        read.xlsx(excel_file, sheet = sheet, rowNames = FALSE, check.names = FALSE)
    })
    names(data_list) <- sheets

    taxonomy_columns <- c("qiime_sklearn", "qiime_confidence", 
                            "vsearch_usearchglobal", "vsearch_identity", 
                            "usearch_utax", "usearch_sintax", "usearch_sintax_80%")
    
    data_fungi <- data_list$fungi
    data_bacteria <- data_list$bacteria
    
    counts_fungi <- data_fungi[, c("OTU_ID", setdiff(names(data_fungi), taxonomy_columns))]
    counts_bacteria <- data_bacteria[, c("OTU_ID", setdiff(names(data_bacteria), taxonomy_columns))]
    
    out$counts <- list(counts_bacteria = counts_bacteria, counts_fungi = counts_fungi)
    
    #-------------------------------------------------------------
    # Tax table
    #-------------------------------------------------------------
    taxa_bacteria <- data_bacteria[, c("OTU_ID", taxonomy_columns)]
    taxa_fungi <- data_fungi[, c("OTU_ID", taxonomy_columns)]
    
    out$tax <- list(taxa_bacteria = taxa_bacteria, taxa_fungi = taxa_fungi)

    #-------------------------------------------------------------
    # Calculate Proportions from Counts
    #-------------------------------------------------------------
    row_sums <- rowSums(counts_bacteria) 
    proportions_bacteria <- sweep(counts_bacteria, 1, row_sums, FUN = "/") 
    proportions_bacteria[is.nan(proportions_bacteria)] <- 0

    row_sums <- rowSums(counts_fungi) 
    proportions_fungi <- sweep(counts_fungi, 1, row_sums, FUN = "/") 
    proportions_fungi[is.nan(proportions_fungi)] <- 0

    out$proportions <- list(proportions_bacteria = proportions_bacteria, proportions_fungi = proportions_fungi)
    }

    # Clean up temporary files
    unlink(temp_dir, recursive = TRUE)


    #-------------------------------------------------------------
    # Read scale information
    #-------------------------------------------------------------
    scale_cols <- c()  
    scale_data <- dataset[, scale_cols]
    out$scale <- scale_data

    #*** Needs to be fixed ***

    #-------------------------------------------------------------
    # Helper function to read files from .zip archives
    #-------------------------------------------------------------

    read_from_zip <- function(zip_path, file_name, is_excel = FALSE, sheet = NULL) {
    temp_dir <- tempdir()
    unzip(zip_path, exdir = temp_dir)
    file_path <- file.path(temp_dir, file_name)
    
    if (is_excel) {
        data <- read.xlsx(file_path, sheet = sheet, check.names = FALSE)
    } else {
        data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    }
    
    unlink(temp_dir, recursive = TRUE) # Clean up
    return(data)
    }

    #-------------------------------------------------------------
    # Read and merge metadata
    #-------------------------------------------------------------

    sex_delivery_zip <- "SI_data3_sex_and_delivery_data.csv.zip"
    diet_data_zip <- "SI_data2_diet_data.csv.zip"
    meds_data_zip <- "SI_data1_allMeds_jan2020.xlsx.zip"

    sex_delivery_data <- read_from_zip(sex_delivery_zip, "SI_data3_sex_and_delivery_data.csv")
    diet_data <- read_from_zip(diet_data_zip, "SI_data2_diet_data.csv")
    meds_data <- read_from_zip(meds_data_zip, "SI_data1_allMeds_jan2020.xlsx", is_excel = TRUE, sheet = 1)

    sex_delivery_data <- sex_delivery_data %>% rename(id = baby_id)

    merged_dataset <- sex_delivery_data %>%
    full_join(diet_data, by = "id") %>%
    full_join(meds_data, by = "id")

    out$metadata <- as.matrix(merged_dataset)

  return(out)
}
