\name{mutt_cache_remove}
\alias{mutt_cache_remove}
\title{Remove saved local study data}
\usage{mutt_cache_remove(studies, cache = NULL)}
\arguments{
\item{studies}{Study IDs to remove, including saved copies from older releases.}
\item{cache}{Directory containing saved local copies (the cache). Defaults to
R's user cache directory for MUTT.}
}
\description{Removes downloaded study copies from your computer. Original input
files downloaded separately with \code{mutt_download_raw()} are kept. A later
call to \code{mutt()} will download the requested study again. Data cannot be
removed while another process is downloading the same study.}
\value{Removed directory paths, invisibly.}
