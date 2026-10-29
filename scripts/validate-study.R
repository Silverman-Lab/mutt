#!/usr/bin/env Rscript
a <- commandArgs(TRUE)
if(length(a)!=1L) stop('Usage: validate-study.R STAGED_STUDY')
source('scripts/release.R');validate_study(a[[1]]);cat('Study validated\n')
