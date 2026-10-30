.mutt_sources <- function() {
  jsonlite::read_json(system.file('extdata','sources.json',package='mutt'))
}

#' Inspect study sources
#'
#' Shows source notes and input information without downloading source files.
#' @param study One study ID.
#' @return Source information, including accession numbers, archived input files,
#'   original notes and known limits on public inputs or processing instructions.
#' @export
mutt_source <- function(study) {
  .mutt_select(study)
  if(length(study)!=1L) stop('Supply one study ID.')
  .mutt_sources()[[study]]
}

#' Find study processing scripts
#'
#' Returns links to the scripts and metadata for the study's recorded source
#' version. These scripts are for maintaining or extending the collection.
#' @param study One study ID.
#' @return Named URLs for study.yml, download.R and process.R.
#' @export
mutt_recipe <- function(study) {
  meta <- mutt_source(study)
  ref <- meta$recipe_commit
  if(is.null(ref) || !grepl('^[a-f0-9]{40}$',ref)) stop('Processing-script links are unavailable: no source version has been recorded.')
  base <- paste0('https://github.com/Silverman-Lab/mutt/blob/',ref,'/studies/',study,'/')
  stats::setNames(paste0(base,c('study.yml','download.R','process.R')),c('metadata','download','process'))
}

#' Download original study inputs
#'
#' Downloads archived inputs, which may include published tables or previously
#' processed data. With reads = TRUE, also downloads public sequencing reads
#' from ENA. These files can be much larger than the processed data returned by
#' [mutt()]. Source notes describe known gaps in public inputs.
#' @param study One study ID.
#' @param dest Destination for source files. Defaults to a sources directory
#'   within MUTT's saved local files (the cache).
#' @param reads Also retrieve public sequencing reads from ENA. Default FALSE.
#' @return Destination path, invisibly. Missing public reads generate warnings.
#' @export
mutt_download_raw <- function(study, dest=NULL, reads=FALSE) {
  if(!is.logical(reads) || length(reads)!=1L || is.na(reads)) stop('reads must be TRUE or FALSE.')
  meta <- mutt_source(study)
  if(is.null(dest)) dest <- file.path(.mutt_cache_root(),'sources',study)
  if(!is.character(dest) || length(dest)!=1L || is.na(dest) || !nzchar(dest)) stop('Invalid destination.')
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  if(!length(meta$archived_inputs)) warning('No archived source inputs are known for ',study)
  for(item in meta$archived_inputs) download_checked(item$url,file.path(dest,item$file),item$sha256)
  if(reads) {
    if(!length(meta$accessions)) warning('No public raw-read accession is known for ',study)
    for(accession in meta$accessions) download_ena_reads(accession,file.path(dest,'reads'))
  }
  invisible(normalizePath(dest))
}
