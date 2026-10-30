\name{mutt_download_raw}
\alias{mutt_download_raw}
\title{Download original study inputs}
\usage{mutt_download_raw(study, dest = NULL, reads = FALSE)}
\arguments{
\item{study}{One study ID.}
\item{dest}{Destination for source files. Defaults to a sources directory within
MUTT's saved local files (the cache).}
\item{reads}{Also retrieve public sequencing reads from ENA. Default FALSE.}
}
\description{Downloads archived inputs, which may include published tables or
previously processed data. With \code{reads = TRUE}, also downloads public
sequencing reads from ENA. These files can be much larger than the processed
data returned by \code{mutt()}. Source notes describe known gaps in public
inputs. Missing public reads generate warnings.}
\value{The destination directory, invisibly.}
