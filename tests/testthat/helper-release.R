fixture_release <- function(ids=c('tiny','second')) {
  root <- tempfile('mutt-fixture-');dir.create(root)
  output <- file.path(root,'assets');dir.create(output)
  rows <- list();paths <- list();objects <- list()
  for(id in ids) {
    staged <- file.path(root,id);dir.create(staged)
    counts <- matrix(c(2,1,3,4),2,dimnames=list(c('S1','S2'),c('ASV1','ASV2')))
    obj <- list(counts=list(original=counts),proportions=list(original=counts/rowSums(counts)),
                tax=list(original=data.frame(Kingdom=c('Bacteria','Bacteria'),row.names=colnames(counts))),
                scale=data.frame(load=c(4,5),row.names=rownames(counts)),
                metadata=data.frame(group=c('a','b'),row.names=rownames(counts)))
    saveRDS(obj,file.path(staged,'study.rds'))
    aligned <- obj;aligned$metadata$group <- c('aligned','aligned')
    saveRDS(aligned,file.path(staged,'aligned.rds'))
    p <- list(study_id=id,schema_version=1L,data_version='1.0.0',selection='fixture',
      study_rds_sha256=mutt:::.mutt_sha256(file.path(staged,'study.rds')),
      variants=list(aligned=list(file='aligned.rds',sha256=mutt:::.mutt_sha256(file.path(staged,'aligned.rds')))))
    jsonlite::write_json(p,file.path(staged,'provenance.json'),auto_unbox=TRUE)
    asset <- file.path(output,paste0(id,'.zip'));old <- setwd(staged)
    status <- utils::zip(asset,c('study.rds','provenance.json','aligned.rds'),flags='-qX',zip=unname(Sys.which('zip')))
    setwd(old);stopifnot(status==0L)
    rows[[id]] <- list(study_id=id,display_name=id,sequencing='amplicon',load='qPCR',status='available',
      release_tag='data-v1.0.0',data_version='1.0.0',
      url=paste0('https://github.com/Silverman-Lab/mutt/releases/download/data-v1.0.0/',id,'.zip'),
      sha256=mutt:::.mutt_sha256(asset),size_bytes=unname(file.info(asset)$size),
      study_rds_sha256=p$study_rds_sha256,provenance_sha256=mutt:::.mutt_sha256(file.path(staged,'provenance.json')))
    paths[[id]] <- asset;objects[[id]] <- obj
  }
  manifest <- file.path(root,'manifest.json')
  jsonlite::write_json(list(schema_version=1L,studies=unname(rows)),manifest,auto_unbox=TRUE)
  list(root=root,manifest=manifest,assets=paths,rows=rows,objects=objects,cache=file.path(root,'cache'))
}
