#!/usr/bin/env Rscript
a <- commandArgs(TRUE)
if(length(a)!=2L) stop('Usage: validate-asset.R ASSET_ZIP STAGED_STUDY')
source('scripts/release.R');validate_asset(a[[1]],a[[2]]);cat('Asset matches staging\n')
