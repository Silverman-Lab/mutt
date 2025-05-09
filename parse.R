# Read metadata
metadata_df <- read_zipped_table(metadata_zip, row.names = NULL)
# Remove empty columns
metadata_df <- metadata_df[, !sapply(metadata_df, function(x) all(is.na(x) | x == ""))]
# Make column names unique
colnames(metadata_df) <- make.unique(colnames(metadata_df))

# Combine Health and status columns
metadata_df$Health.status <- paste(metadata_df$Health, metadata_df$status)
metadata_df <- metadata_df[, !colnames(metadata_df) %in% c("Health", "status")]

# Read SRA data
sra_df <- read_zipped_table(sra_zip, row.names = NULL)
colnames(sra_df) <- c("sampleID", "sra") 