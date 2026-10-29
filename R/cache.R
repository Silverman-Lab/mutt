.mutt_cache_root <- function(cache=NULL) {
  if (is.null(cache)) cache <- getOption('mutt.cache', tools::R_user_dir('mutt','cache'))
  if (!is.character(cache) || length(cache)!=1L || is.na(cache) || !nzchar(cache)) stop('Invalid cache directory.')
  normalizePath(cache,mustWork=FALSE)
}

.mutt_fetch <- function(row, cache, refresh, representation='default') {
  if (row$status != 'available') stop('Study ',row$study_id,' is ',row$status,
    '; no verified public data release is available. ',
    if(!is.null(row$reason)) paste0(row$reason,'. ') else '',
    'See mutt_source() and docs/migration.',call.=FALSE)
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
  if (dir.exists(target) && !file.rename(target,backup)) stop('Could not move the previous cache.')
  if (!file.rename(extracted,target)) {
    if(dir.exists(backup)) file.rename(backup,target)
    stop('Could not atomically publish verified study cache.')
  }
  value
}

.mutt_download_file <- function(url, destination) {
  utils::download.file(url,destination,mode='wb',quiet=TRUE)
}

#' Remove cached processed studies
#' @param studies Character vector of study IDs to remove, across pinned release versions.
#' @param cache Optional cache directory; defaults to the R user cache.
#' @return Removed cache paths, invisibly.
#' @export
mutt_cache_remove <- function(studies, cache=NULL) {
  rows <- .mutt_select(studies)
  root <- .mutt_cache_root(cache)
  paths <- unlist(lapply(names(rows),function(id) {
    candidates <- list.dirs(root,recursive=FALSE,full.names=TRUE)
    unlist(lapply(candidates,function(p) list.files(p,pattern=paste0('^',id,'-[a-f0-9]{64}$'),full.names=TRUE)))
  }),use.names=FALSE)
  for(p in paths) {
    if(dir.exists(paste0(p,'.lock'))) stop('Study cache is in use: ',p)
    unlink(p,recursive=TRUE)
  }
  invisible(paths)
}
