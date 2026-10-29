# Permanent tooling: inputs are canonical staged directories only.
source('R/validation.R',local=TRUE)

validate_study <- function(input) {
  input <- normalizePath(input,mustWork=TRUE)
  value <- .mutt_validate_directory(input)
  p <- jsonlite::read_json(file.path(input,'provenance.json'))
  if(!identical(basename(input),p$study_id)) stop('Staging directory must match provenance study_id.')
  invisible(value)
}

validate_asset <- function(asset, input, row=NULL) {
  expected <- validate_study(input)
  tmp <- tempfile('mutt-verify-');dir.create(tmp)
  on.exit(unlink(tmp,recursive=TRUE),add=TRUE)
  actual <- .mutt_extract_asset(asset,tmp,row)
  for(name in utils::unzip(asset,list=TRUE)$Name) {
    if(!identical(.mutt_sha256(file.path(tmp,name)),.mutt_sha256(file.path(input,name))))
      stop('Release content differs from staged source: ',name)
  }
  if(!identical(actual,expected)) stop('Scientific contents differ from staged source.')
  invisible(TRUE)
}

package_study <- function(input, output, tag) {
  validate_study(input)
  if(!grepl('^data-v[0-9]+[.][0-9]+[.][0-9]+$',tag)) stop('Use a versioned data-vX.Y.Z tag.')
  p <- jsonlite::read_json(file.path(input,'provenance.json'))
  if(!identical(p$data_version,sub('^data-v','',tag))) stop('Staged data version does not match release tag.')
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  archive <- file.path(normalizePath(output),paste0(p$study_id,'.zip'))
  if(file.exists(archive)) stop('Refusing to overwrite asset: ',archive)
  old <- getwd();on.exit(setwd(old),add=TRUE);setwd(input)
  zip_command <- unname(Sys.which('zip'))
  if(!nzchar(zip_command)) stop('The zip command is required for maintainer packaging.')
  members <- unique(c('study.rds','provenance.json',vapply(p$variants,function(v)v$file,character(1))))
  status <- utils::zip(archive,members,flags='-qX9',zip=zip_command)
  if(!identical(status,0L)) stop('ZIP creation failed.')
  setwd(old)
  metadata <- p$metadata
  row <- list(study_id=p$study_id,display_name=metadata$display_name,
    sequencing=metadata$sequencing,load=metadata$load,status='staged',
    release_tag=tag,data_version=p$data_version,schema_version=p$schema_version,
    url=paste0('https://github.com/Silverman-Lab/mutt/releases/download/',tag,'/',p$study_id,'.zip'),
    sha256=.mutt_sha256(archive),size_bytes=unname(file.info(archive)$size),
    study_rds_sha256=.mutt_sha256(file.path(input,'study.rds')),
    provenance_sha256=.mutt_sha256(file.path(input,'provenance.json')))
  validate_asset(archive,input,row)
  row
}

build_release <- function(input, tag, output, manifest=NULL) {
  input <- normalizePath(input,mustWork=TRUE)
  dirs <- sort(list.dirs(input,recursive=FALSE,full.names=TRUE))
  if(!length(dirs)) stop('No staged studies.')
  # Validate the entire input before producing an asset.
  for(d in dirs) {
    validate_study(d)
    invisible(gc(FALSE))
  }
  if(dir.exists(output) && length(list.files(output,all.files=TRUE,no..=TRUE))) stop('Release output must be empty.')
  rows <- lapply(dirs,function(d) package_study(d,output,tag))
  result <- list(schema_version=1L,release_tag=tag,studies=rows)
  jsonlite::write_json(result,file.path(output,'studies.candidate.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
  if(!is.null(manifest)) jsonlite::write_json(result,manifest,pretty=TRUE,auto_unbox=TRUE,null='null')
  invisible(result)
}

.gh <- function(args, output=NULL) {
  if(!nzchar(Sys.which('gh'))) stop('Install the gh CLI and authenticate for maintainer release operations.')
  result <- system2('gh',vapply(args,shQuote,character(1)),stdout=if(is.null(output)) '' else output)
  if(!identical(result,0L)) stop('gh failed; no release verification has been recorded.')
  invisible(result)
}

publish_draft <- function(input, output, repo='Silverman-Lab/mutt') {
  manifest <- jsonlite::read_json(file.path(output,'studies.candidate.json'))
  tag <- manifest$release_tag
  # A successful authenticated listing distinguishes an absent release from an auth/network error.
  listing <- tempfile();on.exit(unlink(listing),add=TRUE)
  .gh(c('release','list','--repo',repo,'--limit','1000','--json','tagName,isDraft'),listing)
  releases <- jsonlite::read_json(listing)
  matching <- Filter(function(x) identical(x$tagName,tag),releases)
  if(length(matching) && !isTRUE(matching[[1]]$isDraft)) stop('Published releases are immutable; choose a new tag.')
  if(!length(matching)) {
    target <- system2('git',c('rev-parse','HEAD'),stdout=TRUE)
    notes <- tempfile();on.exit(unlink(notes),add=TRUE)
    writeLines('Curated processed taxonomic study assets preserving the current mutt API at 01c648b. See per-study provenance and the migration audit for exact default/aligned comparisons and existing parser failures.',notes)
    .gh(c('release','create',tag,'--repo',repo,'--draft','--target',target,
          '--title',tag,'--notes-file',notes))
  }
  # Never replace an asset. If this draft is resumed, verify preexisting assets first.
  existing_path <- tempfile();on.exit(unlink(existing_path),add=TRUE)
  .gh(c('release','view',tag,'--repo',repo,'--json','assets'),existing_path)
  existing <- jsonlite::read_json(existing_path)$assets
  names <- vapply(existing,function(x)x$name,character(1))
  for(row in manifest$studies) {
    asset <- paste0(row$study_id,'.zip');local <- file.path(output,asset)
    validate_asset(local,file.path(input,row$study_id),row)
    if(!asset %in% names) .gh(c('release','upload',tag,local,'--repo',repo))
    retrieved <- tempfile('mutt-download-');dir.create(retrieved)
    tryCatch({
      .gh(c('release','download',tag,'--repo',repo,'--pattern',asset,'--dir',retrieved))
      validate_asset(file.path(retrieved,asset),file.path(input,row$study_id),row)
    },finally=unlink(retrieved,recursive=TRUE))
  }
  jsonlite::write_json(list(release_tag=tag,verified=TRUE,asset_count=length(manifest$studies)),
                       file.path(output,'draft-verification.json'),pretty=TRUE,auto_unbox=TRUE)
  invisible(manifest)
}
