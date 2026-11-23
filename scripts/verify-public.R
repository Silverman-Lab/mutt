#!/usr/bin/env Rscript
# Verify every unauthenticated pinned URL AFTER explicit publication. Does not publish.
a <- commandArgs(TRUE)
if(length(a)!=3L) stop('Usage: verify-public.R STAGING ASSET_DIRECTORY OUTPUT_MANIFEST')
source('scripts/release.R')
manifest <- jsonlite::read_json(file.path(a[[2]],'studies.candidate.json'))
for(i in seq_along(manifest$studies)) {
  row <- manifest$studies[[i]]
  temp <- tempfile(fileext='.zip')
  tryCatch({
    status <- utils::download.file(row$url,temp,mode='wb',quiet=TRUE)
    if(!identical(status,0L)) stop('Public download failed: ',row$study_id)
    validate_asset(temp,file.path(a[[1]],row$study_id),row)
  },finally=unlink(temp))
  manifest$studies[[i]]$status <- 'available'
}
# Retain other curated studies, including explicit unavailable records.
if(file.exists('inst/extdata/studies.json')) {
  existing <- jsonlite::read_json('inst/extdata/studies.json')$studies
  replacement_ids <- vapply(manifest$studies,function(r)r$study_id,character(1))
  retained <- Filter(function(r) !r$study_id %in% replacement_ids,existing)
  manifest$studies <- c(manifest$studies,retained)
}
jsonlite::write_json(manifest,a[[3]],pretty=TRUE,auto_unbox=TRUE,null='null')
