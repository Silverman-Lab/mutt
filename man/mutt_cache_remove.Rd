\name{mutt_cache_remove}
\alias{mutt_cache_remove}
\title{Remove cached processed studies}
\usage{mutt_cache_remove(studies, cache = NULL)}
\arguments{
\item{studies}{Study IDs to remove across cached release versions.}
\item{cache}{Optional cache directory.}
}
\description{Explicitly removes selected processed studies from the R user cache.
Source inputs are retained. Locked caches cannot be removed.}
\value{Removed cache paths, invisibly.}
