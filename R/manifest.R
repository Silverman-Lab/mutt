.mutt_manifest <- function() {
  path <- getOption('mutt.manifest', system.file('extdata', 'studies.json', package='mutt'))
  x <- jsonlite::read_json(path, simplifyVector=FALSE)
  if (!identical(x$schema_version, 1L) || !length(x$studies)) stop('Invalid MUTT manifest.')
  ids <- vapply(x$studies, function(r) r$study_id, character(1))
  if (anyDuplicated(ids) || any(!grepl('^[a-z][a-z0-9_]*$', ids))) stop('Invalid study IDs in manifest.')
  for (r in x$studies) {
    if (!r$status %in% c('available','staged','unavailable')) stop('Invalid release status.')
    if (r$status != 'unavailable') {
      fields <- c('release_tag','url','sha256','size_bytes','data_version',
                  'study_rds_sha256','provenance_sha256')
      if (!all(fields %in% names(r)) || !grepl('^data-v[0-9]+[.][0-9]+[.][0-9]+$',r$release_tag) ||
          !identical(r$url,paste0('https://github.com/Silverman-Lab/mutt/releases/download/',r$release_tag,'/',r$study_id,'.zip')) ||
          !is.numeric(r$size_bytes) || length(r$size_bytes)!=1L || !is.finite(r$size_bytes) || r$size_bytes<=0 ||
          any(!vapply(r[c('sha256','study_rds_sha256','provenance_sha256')],
                      function(s) is.character(s) && length(s)==1L && grepl('^[a-f0-9]{64}$',s),logical(1))))
        stop('Invalid pinned asset for ', r$study_id)
    }
  }
  stats::setNames(x$studies, ids)
}

.mutt_select <- function(studies, manifest=.mutt_manifest()) {
  if (!is.character(studies) || !length(studies) || anyNA(studies) || any(!nzchar(studies)) || anyDuplicated(studies))
    stop('Supply one or more unique study IDs from mutt_studies().', call.=FALSE)
  unknown <- setdiff(studies,names(manifest))
  if (length(unknown)) stop('Unknown study ID(s): ',paste(unknown,collapse=', '),call.=FALSE)
  manifest[studies]
}

#' List curated MUTT studies
#' @return A data frame with identifiers, modalities, release versions and availability.
#' @export
mutt_studies <- function() {
  rows <- .mutt_manifest()
  fields <- c('study_id','display_name','sequencing','load','status','release_tag','data_version','size_bytes')
  do.call(rbind, lapply(rows,function(r) {
    out <- lapply(fields,function(k) if(is.null(r[[k]])) NA else paste(unlist(r[[k]]),collapse='; '))
    names(out) <- fields
    as.data.frame(out,stringsAsFactors=FALSE)
  }))
}
