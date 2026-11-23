#!/usr/bin/env Rscript
# Restore canonical staging from verified, already-built release ZIPs.
# This never calls a study recipe or changes/repackages the assets.
restore_staging <- function(assets, staging) {
  assets <- normalizePath(assets,mustWork=TRUE)
  manifest <- jsonlite::read_json(file.path(assets,'studies.candidate.json'))
  if(!identical(manifest$schema_version,1L) ||
     !grepl('^data-v[0-9]+[.][0-9]+[.][0-9]+$',manifest$release_tag) ||
     !length(manifest$studies)) stop('Invalid candidate release manifest.')
  ids <- vapply(manifest$studies,function(r)r$study_id,character(1))
  if(anyDuplicated(ids) || any(!grepl('^[a-z][a-z0-9_]*$',ids))) stop('Invalid study identifiers.')
  staging <- normalizePath(staging,mustWork=FALSE)
  if(file.exists(staging) || dir.exists(staging)) stop('Refusing to overwrite staging: ',staging)
  dir.create(dirname(staging),recursive=TRUE,showWarnings=FALSE)
  temp <- tempfile('.restore-',tmpdir=dirname(staging));dir.create(temp)
  on.exit(unlink(temp,recursive=TRUE),add=TRUE)
  for(row in manifest$studies) {
    archive <- file.path(assets,paste0(row$study_id,'.zip'))
    if(!file.exists(archive)) stop('Missing asset: ',archive)
    .mutt_extract_asset(archive,file.path(temp,row$study_id),row)
    cat(row$study_id,': verified and restored\n')
    invisible(gc(FALSE))
  }
  if(!file.rename(temp,staging)) stop('Could not publish restored staging directory.')
  invisible(staging)
}

if(sys.nframe()==0L) {
  a <- commandArgs(TRUE)
  if(length(a)!=2L) stop('Usage: restore-staging.R ASSET_DIRECTORY NEW_STAGING_DIRECTORY')
  source('R/validation.R')
  restore_staging(a[[1]],a[[2]])
}
