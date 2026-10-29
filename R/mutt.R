#' Obtain processed microbiome studies
#'
#' Downloads only requested studies from a pinned immutable data release.
#' Every asset is checksum verified and cached in the R user cache directory.
#' All scientific representations inside each stored study are preserved.
#'
#' @param studies One or more study IDs from [mutt_studies()].
#' @param cache Optional cache directory.
#' @param refresh Download and verify again, replacing the cached copy on success.
#' @param representation Stored processing form. Default preserves the old API's
#'   rawdata=FALSE, align_samples=FALSE output. Use "aligned" for its
#'   align_samples=TRUE output, when available.
#' @return A named list of study objects, even for a single study. Each study
#'   contains counts, proportions, tax, scale and metadata; these may contain
#'   named nested representations. No assay branch is selected implicitly.
#' @examples
#' mutt_studies()
#' \dontrun{
#' x <- mutt("vandeputte2017")
#' x <- mutt(c("vandeputte2017", "contijoch2019", "morton2019"))
#' }
#' @export
mutt <- function(studies, cache=NULL, refresh=FALSE, representation='default') {
  if (!is.logical(refresh) || length(refresh)!=1L || is.na(refresh)) stop('refresh must be TRUE or FALSE.')
  if(!is.character(representation) || length(representation)!=1L || is.na(representation) ||
     !grepl('^[a-z][a-z0-9_]*$',representation)) stop('Supply one representation name.')
  rows <- .mutt_select(studies)
  root <- .mutt_cache_root(cache)
  lapply(rows,function(row) .mutt_fetch(row,root,refresh,representation))
}
