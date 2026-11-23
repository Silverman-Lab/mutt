#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if(length(args)!=4L || !identical(args[c(1,3)],c("--input","--output"))) stop("Usage: process.R --input DIRECTORY --output DIRECTORY")
stop("No complete study-specific harmonization recipe existed for this unregistered study. Review archived source workflows before adding a processed release.")
