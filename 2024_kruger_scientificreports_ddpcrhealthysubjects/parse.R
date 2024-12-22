# parse.R

require(readxl)
require(stringr)

parse_2024_kruger_scientificreports_ddpcrhealthysubjects <- function() {
  local <- "2024_kruger_scientificreports_ddpcrhealthysubjects/"  
  
  # Initialize the output list
  out <- list()

  #-------------------------------------------------------------
  # Read counts
  #-------------------------------------------------------------
  counts_file <- paste0(local, "41598_2024_75477_MOESM5_ESM.csv.zip")
  if (file.exists(counts_file)) {
    message("Reading proportions from ", counts_file)
    dataset <- read.table(unz(counts_file, "41598_2024_75477_MOESM5_ESM.csv"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    
    metadata_cols <- c("Subject", "Timepoint", "Milling", "Frequency", "StoolsperDay", 
                    "BristolStoolScale_highest", "WaterContent_perc",
                    "pH", "Calprotectin_ugperg", "MPO_ngperml")
    metadata <- dataset[, metadata_cols]

    counts <- dataset[, c("Subject", "Timepoint", !(colnames(dataset) %in% c(metadata_cols, scale_cols)))]
    out$counts <- counts
    #-------------------------------------------------------------
    # Read scale information
    #-------------------------------------------------------------
    scale_cols <- c("Subject", 
                    "Timepoint",
                    "Mean.Fungi.copies_per.mg.total.weight", 
                    "Mean.Fungi.copies_per.mg.dry.weight", 
                    "Mean.bacteria.copies_per.mg.total.weight", 
                    "Mean.bacteria.copies_per.mg.dry.weight" 
                    )
    scale_data <- dataset[, scale_cols]
    out$scale <- scale_data

    #-------------------------------------------------------------
    # Tax table
    #-------------------------------------------------------------
    taxa <- colnames(dataset[ ,!(colnames(dataset) %in% c(metadata_cols, scale_cols))])

    #-------------------------------------------------------------
    # Calculate Proportions from Counts
    #-------------------------------------------------------------
    row_sums <- rowSums(counts) 
    proportions <- sweep(counts, 1, row_sums, FUN = "/") 
    proportions[is.nan(proportions)] <- 0
    out$proportions <- proportions
  }

  #-------------------------------------------------------------
  # Read and merge metadata
  #-------------------------------------------------------------
  library(tidyverse)
  metadata_file <- paste0(local, "41598_2024_75477_MOESM4_ESM.csv.zip")
  if (file.exists(metadata_file)) {
    message("Reading metadata from ", metadata_file)
    dataset2 <- read.table(unz(metadata_file, "41598_2024_75477_MOESM4_ESM.csv"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    dataset2_t <- dataset2 %>%
    pivot_longer(
        cols = starts_with("Replicate"), 
        names_to = "Replicate",
        values_to = "Value"
    ) %>%
    pivot_wider(          
        names_from = ID,
        values_from = Value
    )
  }
  metadata2_file <- paste0(local, "41598_2024_75477_MOESM3_ESM.csv.zip")
  if (file.exists(metadata2_file)) {
    message("Reading metadata from ", metadata2_file)
    dataset1 <- read.table(unz(metadata2_file, "41598_2024_75477_MOESM3_ESM.csv"), header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
    dataset1 <- dataset1 %>%
    mutate(Replicate = gsub("Replicate ", "Replicate_", ID)) %>%
    select(-ID) 
  }

  merged_dataset <- dataset1 %>% inner_join(dataset2_t, by = "Replicate") 
  merged_dataset <- merged_dataset %>% inner_join(metadata, by = c("Subject","Timepoint"))
  out$metadata <- as.matrix(merged_dataset)

  return(out)
}
