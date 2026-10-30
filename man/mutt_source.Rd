\name{mutt_source}
\alias{mutt_source}
\alias{mutt_recipe}
\title{Inspect study sources and processing scripts}
\usage{
mutt_source(study)
mutt_recipe(study)
}
\arguments{\item{study}{One study ID from \code{mutt_studies()}.}}
\description{Shows source notes and input information without downloading source
files. Processing scripts are for maintaining or extending the collection.}
\value{\code{mutt_source()} returns source information, including accession
numbers, archived input files, original notes and known limits on public inputs
or processing instructions. \code{mutt_recipe()} returns named URLs to
study.yml, download.R and process.R for the study's recorded source version.}
