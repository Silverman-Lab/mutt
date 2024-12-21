# parse.R

require(readxl)
require(stringr)

parse_2022_fromentin_naturemedicine_metacardissubset <- function() {
  local <- "2022_fromentin_naturemedicine_metacardissubset/"  

  out <- list()
  
  #-------------------------------------------------------------
  # Read reduced_feature counts from zipped RData
  #-------------------------------------------------------------
  reduced_feature_zip <- paste0(local, "reduced_feature (1).RData.zip")
  if (file.exists(reduced_feature_zip)) {
    message("Extracting and loading reduced_feature data from ", reduced_feature_zip)
    
    temp_dir <- tempdir()
    unzip(reduced_feature_zip, exdir = temp_dir)
    rdata_files <- list.files(temp_dir, pattern = "\\.RData$", full.names = TRUE)
    
    for (rdata_file in rdata_files) {
      env <- new.env()
      load(rdata_file, envir = env)

      if (exists("reduced_feature", envir = env)) {
        out$counts <- env$reduced_feature
      } else {
        warning("reduced_feature not found in ", rdata_file)
      }
    }
  } else {
    warning("Reduced_feature zip file not found: ", reduced_feature_zip)
  }
  
  #-------------------------------------------------------------
  # Read metamatmetformin metadata from zipped RData
  #-------------------------------------------------------------
  metamatmetformin_zip <- paste0(local, "metaMatMetformin (1).RData.zip")
  if (file.exists(metamatmetformin_zip)) {
    message("Extracting and loading metamatmetformin data from ", metamatmetformin_zip)
    
    temp_dir <- tempdir()
    unzip(metamatmetformin_zip, exdir = temp_dir)
    rdata_files <- list.files(temp_dir, pattern = "\\.RData$", full.names = TRUE)
    
    for (rdata_file in rdata_files) {
      env <- new.env()
      load(rdata_file, envir = env)

      if (exists("metamatmetformin", envir = env)) {
        out$metadata_small <- env$metamatmetformin
      } else {
        warning("metamatmetformin not found in ", rdata_file)
      }
    }
  } else {
    warning("Metamatmetformin zip file not found: ", metamatmetformin_zip)
  }
  
  #-------------------------------------------------------------
  # Read metadata from zipped Excel file
  #-------------------------------------------------------------
  metadata_zip <- paste0(local, "41591_2022_1688_MOESM3_ESM.xlsx.zip")
  if (file.exists(metadata_zip)) {
    message("Extracting and reading metadata from ", metadata_zip)

    temp_dir <- tempdir()
    unzip(metadata_zip, exdir = temp_dir)
    excel_files <- list.files(temp_dir, pattern = "\\.xlsx$", full.names = TRUE)
    
    if (length(excel_files) > 0) {
      metadata_file <- excel_files[1]  
      metadata_sheets <- c("ST10", "ST11", "ST12", "ST13", "ST14")
      
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
  # Checks for reduced_feature and metamatmetformin alignment
  #-------------------------------------------------------------
  if ("counts" %in% names(out) && "metadata_small" %in% names(out)) {
    if (all(rownames(out$counts) %in% rownames(out$metadata_small)) &&
        all(rownames(out$metadata_small) %in% rownames(out$counts))) {
      message("Reduced_feature and metamatmetformin are properly aligned.")
    } else {
      warning("Mismatch detected between reduced_feature rows and metamatmetformin rows.")
    }
  }
  
  #-------------------------------------------------------------
  # Return the combined dataset
  #-------------------------------------------------------------
  return(out)
}