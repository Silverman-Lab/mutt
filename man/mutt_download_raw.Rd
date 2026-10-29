\name{mutt_download_raw}
\alias{mutt_download_raw}
\title{Download original study inputs}
\usage{mutt_download_raw(study, dest = NULL, reads = FALSE)}
\arguments{
\item{study}{One study ID.}
\item{dest}{Destination directory; defaults to a user cache source directory.}
\item{reads}{Also retrieve public sequencing FASTQs from ENA. Default FALSE.}
}
\description{Retrieves exact archived published/preprocessed source inputs from
checksum-pinned historical URLs. Raw reads remain in authoritative repositories
and are downloaded only when explicitly requested. Missing public reads generate
warnings. Source downloads may be large.}
\value{The destination directory, invisibly.}
