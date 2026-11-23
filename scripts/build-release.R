#!/usr/bin/env Rscript
# Usage: Rscript scripts/build-release.R --input .build/data-v1.0.0 --tag data-v1.0.0 --output .build/release-v1 [--draft]
a <- commandArgs(TRUE)
draft <- '--draft' %in% a;a <- a[a!='--draft']
if(length(a)!=6L || !identical(a[c(1,3,5)],c('--input','--tag','--output')))
 stop('Usage: build-release.R --input STAGING --tag data-vX.Y.Z --output ASSETS [--draft]')
source('scripts/release.R')
build_release(a[[2]],a[[4]],a[[6]])
if(draft) publish_draft(a[[2]],a[[6]])
