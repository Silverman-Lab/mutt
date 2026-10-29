\name{mutt_source}
\alias{mutt_source}
\alias{mutt_recipe}
\title{Inspect sources and exact study recipes}
\usage{
mutt_source(study)
mutt_recipe(study)
}
\arguments{\item{study}{One study ID from \code{mutt_studies()}.}}
\description{Inspect original notes, accession numbers, checksums, archived input
URLs, availability limits, and study-specific preprocessing. No source files are
downloaded. Recipe URLs are pinned to a source commit.}
\value{\code{mutt_source()} returns a source metadata list. \code{mutt_recipe()}
returns named URLs to study.yml, download.R and process.R.}
