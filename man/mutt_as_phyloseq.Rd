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
\description{Select one counts table before converting. Metadata must include
every sample ID in the counts; taxonomy must include every feature ID. The
converter matches these tables by their row names and puts them in the counts
table's order. All selected counts are kept. Include microbial-load columns in
metadata if desired. Install the optional phyloseq or TreeSummarizedExperiment
package with BiocManager when using the corresponding converter.}
\value{A phyloseq object, or a TreeSummarizedExperiment usable with mia/miaverse.}
