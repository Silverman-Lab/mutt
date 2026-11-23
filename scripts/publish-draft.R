#!/usr/bin/env Rscript
a <- commandArgs(TRUE)
if(length(a)!=2L) stop('Usage: publish-draft.R STAGING ASSET_DIRECTORY')
source('scripts/release.R');publish_draft(a[[1]],a[[2]])
