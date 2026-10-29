.mutt_conversion_tables <- function(counts, metadata, taxonomy) {
  if(!is.matrix(counts) && !is.data.frame(counts)) stop('Select a single counts table from the study branches.')
  counts <- as.matrix(counts)
  if(!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) || any(counts<0)) stop('Counts must be finite nonnegative numeric values.')
  for(ids in dimnames(counts)) if(is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop('Counts need unique sample and feature identifiers.')
  if(is.null(rownames(counts)) || is.null(colnames(counts))) stop('Counts need sample and feature identifiers.')
  align_table <- function(x,ids,label) {
    if(is.null(x)) return(data.frame(row.names=ids))
    x <- as.data.frame(x)
    if(anyDuplicated(rownames(x)) || !all(ids %in% rownames(x))) stop(label,' row names must include all selected identifiers. Set them explicitly; no samples or features are dropped.')
    x[ids,,drop=FALSE]
  }
  list(counts=counts,metadata=align_table(metadata,rownames(counts),'Metadata'),
       taxonomy=align_table(taxonomy,colnames(counts),'Taxonomy'))
}

#' Convert a selected counts table to phyloseq or miaverse
#'
#' Select a scientific branch explicitly before converting. Counts have samples
#' in rows and features in columns. Set metadata/taxonomy row names to their
#' identifiers; all counts are retained and tables are reordered by identifier.
#' Include microbial-load columns in metadata if desired. Conversion never
#' substitutes proportions for counts or infers missing taxonomy.
#' @param counts One numeric count matrix or data frame.
#' @param metadata Optional sample table, with sample IDs as row names.
#' @param taxonomy Optional feature table, with feature IDs as row names.
#' @return A phyloseq object or a TreeSummarizedExperiment for use with mia/miaverse.
#' @export
mutt_as_phyloseq <- function(counts, metadata=NULL, taxonomy=NULL) {
  x <- .mutt_conversion_tables(counts,metadata,taxonomy)
  if(!requireNamespace('phyloseq',quietly=TRUE)) stop('Install phyloseq with BiocManager::install("phyloseq") to convert.')
  parts <- list(phyloseq::otu_table(x$counts,taxa_are_rows=FALSE))
  if(ncol(x$metadata)) parts <- c(parts,list(phyloseq::sample_data(x$metadata)))
  if(ncol(x$taxonomy)) parts <- c(parts,list(phyloseq::tax_table(as.matrix(x$taxonomy))))
  do.call(phyloseq::phyloseq,parts)
}

#' @rdname mutt_as_phyloseq
#' @export
mutt_as_miaverse <- function(counts, metadata=NULL, taxonomy=NULL) {
  x <- .mutt_conversion_tables(counts,metadata,taxonomy)
  if(!requireNamespace('TreeSummarizedExperiment',quietly=TRUE) || !requireNamespace('S4Vectors',quietly=TRUE))
    stop('Install TreeSummarizedExperiment with BiocManager::install("TreeSummarizedExperiment") to convert.')
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays=list(counts=t(x$counts)),colData=S4Vectors::DataFrame(x$metadata),
    rowData=S4Vectors::DataFrame(x$taxonomy))
}
