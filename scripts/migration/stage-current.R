#!/usr/bin/env Rscript
# ONE-TIME ONLY: copy captured current-API objects into permanent staging.
a <- commandArgs(TRUE)
if(length(a)!=3L) stop('Usage: stage-current.R DEFAULT_OBJECTS ALIGNED_OBJECTS STAGING')
metadata <- jsonlite::read_json('inst/extdata/sources.json')
artifacts <- read.csv('docs/migration/baseline-artifacts.csv',stringsAsFactors=FALSE)
report <- list()
for(id in names(metadata)) {
  meta <- metadata[[id]];legacy <- meta$legacy_id
  original <- file.path(a[[1]],paste0(legacy,'.rds'))
  if(!file.exists(original)) next
  dest <- file.path(a[[3]],id)
  if(dir.exists(dest)) stop('Refusing to overwrite staged study ',id)
  dir.create(dest,recursive=TRUE)
  stopifnot(file.copy(original,file.path(dest,'study.rds')))
  sha <- function(p) digest::digest(p,algo='sha256',file=TRUE)
  stopifnot(identical(sha(original),sha(file.path(dest,'study.rds'))))
  value <- readRDS(original)
  stopifnot(identical(value,readRDS(file.path(dest,'study.rds'))))
  variants <- list()
  aligned <- file.path(a[[2]],paste0(legacy,'.rds'))
  if(file.exists(aligned)) {
    stopifnot(file.copy(aligned,file.path(dest,'aligned.rds')))
    stopifnot(identical(sha(aligned),sha(file.path(dest,'aligned.rds'))))
    variants$aligned <- list(file='aligned.rds',sha256=sha(aligned),
      legacy_call='mutt(study, rawdata=FALSE, align_samples=TRUE, functional=FALSE)')
  }
  p <- list(study_id=id,schema_version=1L,data_version='1.0.0',
    baseline_status='current API output preserved exactly, excluding removed functional field',
    legacy_commit='01c648b1f736b1af019d1c0f0231b73de2c21e49',legacy_path=legacy,
    legacy_object_call='mutt(study, rawdata=FALSE, align_samples=FALSE, functional=FALSE)[[1]]',
    captured_object_sha256=sha(original),study_rds_sha256=sha(original),
    selection='All branches from current mutt(); no count, metadata, taxon-name, load, phylogeny or demographic changes.',
    conversion='Remove only the empty inferred-function field; serialize once, then byte-copy the captured RDS.',
    variants=variants,metadata=meta[c('display_name','sequencing','load','pmid','accessions')],
    legacy_artifacts=artifacts[artifacts$study==legacy,c('path','sha256','size_bytes')])
  jsonlite::write_json(p,file.path(dest,'provenance.json'),pretty=TRUE,auto_unbox=TRUE,null='null')
  report[[id]] <- data.frame(study_id=id,legacy_id=legacy,default_identical=TRUE,
    aligned_identical=length(variants)>0,default_sha256=sha(original),
    scientific_sha256=digest::digest(value,algo='sha256'))
  cat(id,': exact copy\n')
}
write.csv(do.call(rbind,report),'docs/migration/current-object-comparison.csv',row.names=FALSE)
