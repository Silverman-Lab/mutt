#' Obtain processed microbiome studies
#'
#' Downloads only the processed studies you request and saves a local copy
#' for reuse. Only studies marked "available" by [mutt_studies()] can be
#' downloaded. You do not need to run the study's processing scripts.
#'
#' @param studies One or more study IDs from [mutt_studies()].
#' @param cache Directory for saved local copies (the cache). If NULL, uses
#'   `getOption("mutt.cache")` or, if unset, R's user cache directory for \emph{mutt}.
#' @param refresh Download and check again. A failed download leaves the previous
#'   saved copy usable.
#' @param representation Which saved version to load. "default" keeps the
#'   original sample groups, which can differ between tables. "aligned" loads
#'   a second saved version, when provided, prepared by matching sample IDs.
#'   This is sample matching, not DNA sequence alignment.
#' @details For example, if counts contain samples A, B and C but load measurements
#'   exist only for A and C, the "aligned" version may limit counts to A and C.
#'   Rules vary by study; tables need not contain the same samples or row order.
#'   Match sample IDs when combining tables, rather than using row positions.
#'   Both choices load saved tables; [mutt()] does no further sample filtering.
#'
#'   Each download is checked automatically against its expected file size and
#'   SHA-256 checksum, a fingerprint of the file contents. Later requests reuse
#'   the checked local copy.
#' @return A named list of study objects, even for a single study. Each study
#'   contains `counts` (sequence counts), `proportions` (relative abundances),
#'   `tax` (taxonomic annotations), `scale` (microbial-load measurements with
#'   study-specific names and units) and `metadata` (sample information). A field
#'   may contain several named tables for different assays, classifiers or
#'   cohorts. Select the table you need for your analysis.
#' @examples
#' mutt_studies()
#' \dontrun{
#' x <- mutt("vandeputte2017")
#' x <- mutt(c("vandeputte2017", "contijoch2019", "morton2019"))
#' }
#' @seealso [mutt_studies()], [mutt_source()], [mutt_cache_remove()]
#' @export
mutt <- function(studies, cache=NULL, refresh=FALSE, representation='default') {
  if (!is.logical(refresh) || length(refresh)!=1L || is.na(refresh)) stop('refresh must be TRUE or FALSE.')
  if(!is.character(representation) || length(representation)!=1L || is.na(representation) ||
     !grepl('^[a-z][a-z0-9_]*$',representation)) stop('Supply one saved-version name, such as default or aligned.')
  rows <- .mutt_select(studies)
  root <- .mutt_cache_root(cache)
  lapply(rows,function(row) .mutt_fetch(row,root,refresh,representation))
}
