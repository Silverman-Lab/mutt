#!/usr/bin/env Rscript
# ONE-TIME ONLY. Run with the PRE-migration mutt installed in an isolated library.
# Arguments: legacy library, legacy studies directory, parser index, output directory.
a <- commandArgs(TRUE)
if(length(a)!=4L) stop('Usage: capture-current.R OLD_R_LIBRARY LEGACY_STUDIES PARSER_INDEX OUTPUT')
library(mutt,lib.loc=normalizePath(a[[1]],mustWork=TRUE))
if(!'functional' %in% names(formals(mutt))) stop('Expected the pre-migration mutt API.')
if(requireNamespace('dada2',quietly=TRUE)) for(fn in c('assignTaxonomy','assignSpecies','dada','filterAndTrim'))
  trace(fn,tracer=quote(stop('Migration refuses scientific reprocessing')),where=asNamespace('dada2'),print=FALSE)
ids <- read.csv(a[[3]],stringsAsFactors=FALSE)$study
output <- normalizePath(a[[4]],mustWork=FALSE)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
report <- list()
for(aligned in c(FALSE,TRUE)) for(id in ids) {
  mode <- if(aligned) 'aligned' else 'default'
  dest <- file.path(output,mode,paste0(id,'.rds'))
  if(file.exists(dest)) stop('Refusing to overwrite baseline object: ',dest)
  result <- suppressWarnings(mutt(id,base_directory=normalizePath(a[[2]],mustWork=TRUE),
                                  rawdata=FALSE,align_samples=aligned,functional=FALSE))
  audit <- attr(result,'audit');audit$representation <- mode
  if(identical(audit$status[[1]],'success')) {
    value <- result[[1]];value[['function']] <- NULL # intentionally removed subsystem only
    dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
    saveRDS(value,dest,version=3)
    stopifnot(identical(value,readRDS(dest)))
  }
  report[[paste(id,mode)]] <- audit
  write.csv(do.call(rbind,report),file.path(output,'capture-report.csv'),row.names=FALSE)
}
