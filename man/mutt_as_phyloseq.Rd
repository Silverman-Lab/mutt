\name{mutt_as_phyloseq}
\alias{mutt_as_phyloseq}
\alias{mutt_as_miaverse}
\title{Convert a selected counts table to phyloseq or miaverse}
\usage{
mutt_as_phyloseq(counts, metadata = NULL, taxonomy = NULL)
mutt_as_miaverse(counts, metadata = NULL, taxonomy = NULL)
}
\arguments{
\item{counts}{One numeric count matrix or data frame, with samples as rows and
features as columns. Both dimensions require unique nonempty identifiers.}
\item{metadata}{Optional sample table with sample IDs as row names. Include load
measurements as columns here to carry them into the converted object.}
\item{taxonomy}{Optional feature table with feature IDs as row names.}
}
\description{Choose a scientific branch explicitly. Conversion retains all counts
and aligns tables by identifiers without selecting taxa or samples implicitly.
Install the optional Bioconductor phyloseq or TreeSummarizedExperiment package
only when using the corresponding conversion.}
\value{A phyloseq object, or a TreeSummarizedExperiment usable with mia/miaverse.}
