.mutt_cache_root <- function(cache=NULL) {
  if (is.null(cache)) cache <- getOption('mutt.cache', tools::R_user_dir('mutt','cache'))
  if (!is.character(cache) || length(cache)!=1L || is.na(cache) || !nzchar(cache)) stop('Invalid directory for saved study copies.')
  normalizePath(cache,mustWork=FALSE)
}

.mutt_fetch <- function(row, cache, refresh, representation='default') {
  if (row$status != 'available') stop('Study ',row$study_id,' is ',row$status,
    '; this study is not yet available for download. ',
    if(!is.null(row$reason)) paste0(row$reason,' ') else '',
    'See mutt_source("',row$study_id,'") for sources and availability.',call.=FALSE)
  parent <- file.path(cache,row$release_tag)
  target <- file.path(parent,paste0(row$study_id,'-',row$sha256))
  if (dir.exists(target) && !refresh) return(.mutt_validate_directory(target,row,representation))
  dir.create(parent,recursive=TRUE,showWarnings=FALSE)
  lock <- paste0(target,'.lock')
  if (!dir.create(lock,showWarnings=FALSE)) stop('Another process is preparing ',row$study_id,'. Retry after it finishes.')
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  stage <- tempfile('.download-',tmpdir=parent);dir.create(stage)
  on.exit(unlink(stage,recursive=TRUE),add=TRUE)
  archive <- file.path(stage,'asset.zip')
  status <- .mutt_download_file(row$url,archive)
  if (!identical(status,0L)) stop('Study download failed.')
  extracted <- file.path(stage,'verified')
  .mutt_extract_asset(archive,extracted,row)
  value <- .mutt_validate_directory(extracted,row,representation)
  # Keep a verified existing cache usable if refresh fails. Both renames share a filesystem.
  backup <- file.path(stage,'previous')
  if (dir.exists(target) && !file.rename(target,backup)) stop('Could not move the previous saved study copy.')
  if (!file.rename(extracted,target)) {
    if(dir.exists(backup)) file.rename(backup,target)
    stop('Could not save the verified study copy.')
  }
  value
}

.mutt_download_file <- function(url, destination) {
  utils::download.file(url,destination,mode='wb',quiet=TRUE)
}

#' Remove saved local study data
#'
#' Removes downloaded study copies from your computer. Original input files
#' downloaded separately with [mutt_download_raw()] are kept. A later call to
#' [mutt()] will download the requested study again. Removal stops if the saved
#' copy is in use by another download.
#' @param studies Study IDs to remove, including saved copies from older releases.
#' @param cache Directory containing saved local copies (the cache). If NULL,
#'   uses `getOption("mutt.cache")` or, if unset, R's user cache directory for \emph{mutt}.
#' @return Removed directory paths, invisibly.
#' @export
mutt_cache_remove <- function(studies, cache=NULL) {
  rows <- .mutt_select(studies)
  root <- .mutt_cache_root(cache)
  paths <- unlist(lapply(names(rows),function(id) {
    candidates <- list.dirs(root,recursive=FALSE,full.names=TRUE)
    unlist(lapply(candidates,function(p) list.files(p,pattern=paste0('^',id,'-[a-f0-9]{64}$'),full.names=TRUE)))
  }),use.names=FALSE)
  for(p in paths) {
    if(dir.exists(paste0(p,'.lock'))) stop('Saved study copy is in use: ',p)
    unlink(p,recursive=TRUE)
  }
  invisible(paths)
}
