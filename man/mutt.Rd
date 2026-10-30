\name{mutt}
\alias{mutt}
\title{Obtain processed microbiome studies}
\usage{mutt(studies, cache = NULL, refresh = FALSE, representation = "default")}
\arguments{
\item{studies}{One or more study IDs from \code{mutt_studies()}.}
\item{cache}{Directory for saved local copies (the cache). Defaults to R's user
cache directory for MUTT.}
\item{refresh}{Download and check again. A failed download leaves the previous
saved copy usable.}
\item{representation}{Which saved version to load. \code{"default"} keeps the
original sample groups, which can differ between tables. \code{"aligned"} loads
a second saved version using the previous package's study-specific rules for
matching sample IDs. This is sample matching, not DNA sequence alignment.}
}
\description{Downloads only the processed studies you request and saves a local
copy for reuse. Only studies marked \code{available} by \code{mutt_studies()}
can be downloaded. You do not need to run the study's processing scripts.}
\details{For example, if counts contain samples A, B and C but load measurements
exist only for A and C, the \code{"aligned"} version may limit counts to A and C.
Rules vary by study; tables need not contain the same samples or row order.
Match sample IDs when combining tables, rather than using row positions.

Each download is checked automatically against its expected file size and
SHA-256 checksum, a fingerprint of the file contents. Later requests reuse the
checked local copy.}
\value{A named list of study objects, even for a single study. Each study contains
counts, proportions, tax, scale and metadata; these may contain several named
tables for different assays, classifiers or cohorts. Select the table you need
for your analysis.}
\examples{
mutt_studies()
\dontrun{
x <- mutt("vandeputte2017")
x <- mutt(c("vandeputte2017", "contijoch2019", "morton2019"))
}
}
\seealso{\code{\link{mutt_source}}, \code{\link{mutt_cache_remove}}}
