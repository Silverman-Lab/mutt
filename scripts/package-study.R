#!/usr/bin/env Rscript
a <- commandArgs(TRUE)
if(length(a)!=3L) stop('Usage: package-study.R STAGED_STUDY OUTPUT_DIRECTORY data-vX.Y.Z')
source('scripts/release.R');print(package_study(a[[1]],a[[2]],a[[3]]))
