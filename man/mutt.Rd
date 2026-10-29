\name{mutt}
\alias{mutt}
\title{Obtain processed microbiome studies}
\usage{mutt(studies, cache = NULL, refresh = FALSE, representation = "default")}
\arguments{
\item{studies}{One or more study IDs from \code{mutt_studies()}.}
\item{cache}{Optional cache directory, defaulting to the R user cache.}
\item{refresh}{Download and verify again. Failed refresh preserves the previous copy.}
\item{representation}{Stored processing form: \code{"default"} preserves the old
rawdata=FALSE, align_samples=FALSE result. \code{"aligned"} preserves its
align_samples=TRUE result, when available.}
}
\description{Downloads only requested datasets from pinned versioned GitHub
Release assets. Checks SHA-256 and byte size, verifies extraction, and caches
completed datasets in the R user cache. Staged assets are not public downloads.}
\value{A named list of scientific study objects, even for one study. Objects retain
counts, proportions, tax, scale, metadata and heterogeneous nested branches.}
\examples{
mutt_studies()
\dontrun{x <- mutt("vandeputte2017")}
}
\seealso{\code{\link{mutt_source}}, \code{\link{mutt_cache_remove}}}
