\name{mutt_studies}
\alias{mutt_studies}
\title{List curated microbiome studies}
\usage{mutt_studies()}
\description{Lists the study IDs, sequencing and load-measurement methods,
data versions and availability recorded in the package's study catalog. It does
not download study data. Status \code{available} means a verified public download
can be requested; \code{staged} means data are prepared for release; and
\code{unavailable} means this release has no usable processed object.}
\value{A data frame containing study_id, display_name, sequencing, load, status,
release_tag, data_version and size_bytes.}
