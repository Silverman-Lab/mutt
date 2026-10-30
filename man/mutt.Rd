\name{mutt}
\alias{mutt}
\title{Obtain processed microbiome studies}
\usage{mutt(studies, cache = NULL, refresh = FALSE, representation = "default")}
\arguments{
\item{studies}{One or more study IDs from \code{mutt_studies()}.}
\item{cache}{Optional cache directory, defaulting to the R user cache.}
\item{refresh}{Download and verify again. Failed refresh preserves the previous copy.}
\item{representation}{Which stored form to load. \code{"default"} retains the
original sample sets, which can differ between tables. \code{"aligned"} retrieves
the previous package's study-specific sample-aligned form, when available.}
}
\description{Retrieves saved study objects from a specific versioned GitHub
Release. Only requested studies are downloaded. Each download is checked against
its recorded file size and SHA-256 checksum, a fingerprint of the file's exact
contents. Checked data are saved in R's user cache for reuse. A study marked
\code{staged} in the catalog is prepared for release but cannot yet be downloaded
through this function.}
\value{A named list of scientific study objects, even for one study. Objects retain
counts, proportions, tax, scale, metadata and heterogeneous nested branches.}
\examples{
mutt_studies()
\dontrun{x <- mutt("vandeputte2017")}
}
\seealso{\code{\link{mutt_source}}, \code{\link{mutt_cache_remove}}}
