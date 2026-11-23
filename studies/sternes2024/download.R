#!/usr/bin/env Rscript
# Study-specific source inputs and accession choices are in this directory's study.yml.
# Default: retrieve the exact stored inputs, including published/preprocessed tables.
# --reads additionally retrieves public FASTQs from ENA; it may be very large.
argv <- commandArgs(trailingOnly=TRUE)
if (length(argv) < 2L || argv[[1]] != "--dest" ||
    !all(argv[-c(1L,2L)] %in% "--reads"))
  stop("Usage: Rscript download.R --dest /source-directory [--reads]")
script <- normalizePath(sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[[1]]))
source(file.path(dirname(script),"..","..","scripts","download-helpers.R"))
meta <- jsonlite::read_json(file.path(dirname(script),"study.yml"))
dest <- argv[[2]]
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
for (item in meta$archived_inputs) download_checked(item$url, file.path(dest,item$file), item$sha256)
if ("--reads" %in% argv) {
  if (!length(meta$accessions)) warning("No public raw-read accession is known for this study.")
  for (accession in meta$accessions) download_ena_reads(accession, file.path(dest,"reads"))
}
message("Source availability and original preprocessing notes: ",file.path(dirname(script),"study.yml"))
