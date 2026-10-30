\name{mutt_studies}
\alias{mutt_studies}
\title{List curated MUTT studies}
\usage{mutt_studies()}
\description{Lists studies without downloading their data. Status
\code{available} means processed data can be downloaded; \code{staged} means
the data are prepared but cannot yet be downloaded; and \code{unavailable}
means this release provides no processed data for the study.}
\value{A data frame containing study_id, display_name, sequencing, load, status,
release_tag, data_version and size_bytes. These describe study IDs, names,
sequencing and load-measurement methods, availability, data versions, release
tags and download sizes.}
