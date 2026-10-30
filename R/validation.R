.mutt_sha256 <- function(path) digest::digest(path, algo='sha256', file=TRUE)

# Preserve heterogeneous study branches. Validation does not normalize scientific values.
.mutt_validate_object <- function(x) {
  required <- c('counts','proportions','tax','scale','metadata')
  if (!is.list(x) || is.data.frame(x) || is.null(names(x)) || anyDuplicated(names(x)) ||
      !all(required %in% names(x))) stop('Invalid study object: expected a named list containing counts, proportions, tax, scale and metadata.')
  if ('function' %in% names(x)) stop('Inferred functional profiles are outside the MUTT schema.')
  invisible(x)
}

.mutt_validate_directory <- function(path, row=NULL, representation='default') {
  files <- file.path(path,c('study.rds','provenance.json'))
  if (!all(file.exists(files))) stop('Incomplete study directory.')
  p <- jsonlite::read_json(files[[2]])
  if (!identical(p$schema_version,1L) || is.null(p$study_id) ||
      !grepl('^[a-z][a-z0-9_]*$',p$study_id) ||
      is.null(p$data_version) || !nzchar(p$data_version) ||
      is.null(p$selection) || !nzchar(p$selection) ||
      !identical(p$study_rds_sha256,.mutt_sha256(files[[1]]))) stop('Invalid study source record or data checksum.')
  if (!is.null(row)) {
    if (!identical(p$study_id,row$study_id) || !identical(p$data_version,row$data_version) ||
        !identical(.mutt_sha256(files[[1]]),row$study_rds_sha256) ||
        !identical(.mutt_sha256(files[[2]]),row$provenance_sha256)) stop('Saved study copy does not match the version in the study catalog.')
  }
  variants <- p$variants
  if(length(variants)) for(name in names(variants)) {
    v <- variants[[name]]
    if(!grepl('^[a-z][a-z0-9_]*$',name) || !grepl('^[a-z][a-z0-9_]*[.]rds$',v$file) ||
       !identical(.mutt_sha256(file.path(path,v$file)),v$sha256)) stop('Invalid checksum for a saved study form.')
    .mutt_validate_object(readRDS(file.path(path,v$file)))
  }
  selected <- files[[1]]
  if(representation != 'default') {
    if(!representation %in% names(variants)) stop('Requested study form is unavailable: ',representation)
    selected <- file.path(path,variants[[representation]]$file)
  }
  x <- readRDS(selected)
  .mutt_validate_object(x)
  x
}

.mutt_extract_asset <- function(archive, destination, row=NULL) {
  if (!is.null(row) && (!identical(unname(file.info(archive)$size),as.numeric(row$size_bytes)) ||
      !identical(.mutt_sha256(archive),row$sha256))) stop('Downloaded file size or SHA-256 verification failed.')
  entries <- utils::unzip(archive,list=TRUE)$Name
  if (!all(c('study.rds','provenance.json') %in% entries) || anyDuplicated(entries) ||
      any(!grepl('^(provenance[.]json|[a-z][a-z0-9_]*[.]rds)$',entries)))
    stop('Invalid archive: required files are missing or file names are unexpected or repeated.')
  dir.create(destination,recursive=TRUE,showWarnings=FALSE)
  utils::unzip(archive,exdir=destination)
  value <- .mutt_validate_directory(destination,row)
  p <- jsonlite::read_json(file.path(destination,'provenance.json'))
  expected <- unique(c('study.rds','provenance.json',vapply(p$variants,function(v)v$file,character(1))))
  if(!setequal(entries,expected)) stop('Invalid archive: files do not match the study source record.')
  value
}
