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
jsonlite::write_json(manifest,a[[3]],pretty=TRUE,auto_unbox=TRUE,null='null')
