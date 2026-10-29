#!/usr/bin/env Rscript
# ONE-TIME ONLY: extract intact scientific objects from the historical saved library.
# Run from repository root. This never calls a parser, DADA2 or MetaPhlAn.
args <- commandArgs(TRUE)
source_file <- if(length(args)) args[[1]] else '.build/legacy/datasetsfromrepo.RData'
staging <- if(length(args)>1) args[[2]] else '.build/data-v1.0.0'
e <- new.env(parent=emptyenv());load(source_file,envir=e)
x <- e$microbialscalerepository
if(!is.list(x)) stop('Expected microbialscalerepository in historical artifact')
map <- read.csv('scripts/migration/v1-selection.csv',stringsAsFactors=FALSE)
legacy_sha <- digest::digest(source_file,algo='sha256',file=TRUE)
for(i in seq_len(nrow(map))) {
 r <- map[i,]; value <- x[[r$output_name]]
 if(is.null(value)) stop('Missing historical object: ',r$output_name)
 dest <- file.path(staging,r$study_id)
 if(dir.exists(dest)) stop('Refusing to overwrite staged study: ',dest)
 dir.create(dest,recursive=TRUE)
 saveRDS(value,file.path(dest,'study.rds'),version=3)
 restored <- readRDS(file.path(dest,'study.rds'))
 if(!identical(value,restored)) stop('Scientific content changed for ',r$study_id)
 metadata <- jsonlite::read_json(file.path('studies',r$study_id,'study.yml'),simplifyVector=FALSE)
 provenance <- list(study_id=r$study_id,schema_version=1L,data_version='1.0.0',
  baseline_status='provisional: manuscript not supplied',
  legacy_commit='53a5406',legacy_git_blob='5ef175e49c7eec8b4725cf33ac5aa2bf087bbbc8',
  legacy_path='datasetsfromrepo.RData',legacy_sha256=legacy_sha,
  legacy_object=paste0('microbialscalerepository$',r$output_name),
  selection='Intact stored study object; preserve all original, reprocessed, classifier, cohort and load branches. No preprocessing rerun.',
  conversion='Extract list entry and serialize as RDS; identical() after round trip.',
  scientific_sha256=digest::digest(value,algo='sha256'),
  study_rds_sha256=digest::digest(file.path(dest,'study.rds'),algo='sha256',file=TRUE),
  metadata=metadata[c('display_name','sequencing','load','pmid','accessions')])
 jsonlite::write_json(provenance,file.path(dest,'provenance.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
 cat(r$study_id,': identical scientific object\n')
}
