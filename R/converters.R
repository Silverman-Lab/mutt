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
#' Select one counts table before converting. Counts have samples in rows and
#' features in columns. Metadata must include every sample ID in the counts;
#' taxonomy must include every feature ID. The converter matches these tables
#' by their row names and puts them in the counts table's order. All selected
#' counts are kept. Include microbial-load columns in metadata if desired.
#' Install the optional phyloseq or TreeSummarizedExperiment package with
#' BiocManager when using the corresponding converter.
#' @param counts One numeric count matrix or data frame, with unique sample IDs
#'   as row names and unique feature IDs as column names.
#' @param metadata Optional sample table, with sample IDs as row names.
#' @param taxonomy Optional feature table, with feature IDs as row names.
#' @return A phyloseq object or a TreeSummarizedExperiment for use with mia/miaverse.
#' @export
mutt_as_phyloseq <- function(counts, metadata=NULL, taxonomy=NULL) {
  x <- .mutt_conversion_tables(counts,metadata,taxonomy)
  if(!requireNamespace('phyloseq',quietly=TRUE)) stop('Install phyloseq with BiocManager::install("phyloseq") to convert.')
  # A single OTU component is returned unchanged by phyloseq's constructor.
  # Supply the existing sample identifiers when no covariates were provided.
  if(!ncol(x$metadata)) x$metadata <- data.frame(sample_id=rownames(x$counts),row.names=rownames(x$counts))
  parts <- list(phyloseq::otu_table(x$counts,taxa_are_rows=FALSE),phyloseq::sample_data(x$metadata))
  if(ncol(x$taxonomy)) parts <- c(parts,list(phyloseq::tax_table(as.matrix(x$taxonomy))))
  result <- do.call(phyloseq::phyloseq,parts)
  if(!identical(phyloseq::sample_names(result),rownames(x$counts)) ||
     !identical(phyloseq::taxa_names(result),colnames(x$counts)))
    stop('phyloseq changed sample or feature identifiers. Use miaverse or provide compatible identifiers explicitly.')
  result
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
