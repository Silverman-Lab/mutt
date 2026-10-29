# Trivial retrieval helpers only; study choices belong in study.yml/download.R.
download_checked <- function(url, path, sha256) {
  if (file.exists(path) && identical(digest::digest(path,algo='sha256',file=TRUE),sha256)) return(invisible(path))
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  temp <- tempfile('.source-',tmpdir=dirname(path));on.exit(unlink(temp),add=TRUE)
  status <- utils::download.file(url,temp,mode='wb',quiet=TRUE)
  if (!identical(status,0L) || !identical(digest::digest(temp,algo='sha256',file=TRUE),sha256))
    stop('Source download or checksum verification failed: ',url)
  if(file.exists(path)) stop('Refusing to overwrite a different source file: ',path)
  if(!file.rename(temp,path)) stop('Could not publish source file: ',path)
  invisible(path)
}

download_ena_reads <- function(accession, dest) {
  url <- paste0('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=',accession,
    '&result=read_run&fields=run_accession,fastq_ftp,fastq_md5&format=tsv')
  runs <- tryCatch(utils::read.delim(url,stringsAsFactors=FALSE),error=function(e) NULL)
  if(is.null(runs) || !nrow(runs) || !all(c('fastq_ftp','fastq_md5') %in% names(runs))) {
    warning('No publicly downloadable ENA FASTQs found for ',accession);return(invisible(NULL))
  }
  dir.create(dest,recursive=TRUE,showWarnings=FALSE)
  for(i in seq_len(nrow(runs))) {
    paths <- strsplit(runs$fastq_ftp[[i]],';',fixed=TRUE)[[1]]
    hashes <- strsplit(runs$fastq_md5[[i]],';',fixed=TRUE)[[1]]
    if(!length(paths) || length(paths)!=length(hashes)) {warning('Missing FASTQs/checksums for ',runs$run_accession[[i]]);next}
    for(j in seq_along(paths)) {
      path <- file.path(dest,basename(paths[[j]]));tmp <- tempfile('.fastq-',tmpdir=dest)
      tryCatch({
        if(file.exists(path) && identical(unname(tools::md5sum(path)),hashes[[j]])) next
        if(file.exists(path)) stop('Existing FASTQ checksum mismatch: ',path)
        status <- utils::download.file(paste0('https://',sub('^ftp://','',paths[[j]])),tmp,mode='wb')
        if(!identical(status,0L) || !identical(unname(tools::md5sum(tmp)),hashes[[j]])) stop('FASTQ checksum mismatch')
        if(!file.rename(tmp,path)) stop('Could not publish FASTQ')
      },finally=unlink(tmp))
    }
  }
  invisible(dest)
}
