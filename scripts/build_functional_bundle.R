#!/usr/bin/env Rscript

options <- list(
  optparse::make_option(
    "--study-dir",
    type = "character",
    help = "Study directory containing a validated functional/ cache."
  ),
  optparse::make_option(
    "--max-asset-bytes",
    type = "double",
    default = 1900000000,
    help = "Maximum size of any published file [default: %default]."
  ),
  optparse::make_option(
    "--target-shard-bytes",
    type = "double",
    default = 268435456,
    help = "Target uncompressed contribution-shard size [default: %default]."
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = options))
if (is.null(args$study_dir) || !nzchar(args$study_dir)) {
  stop("--study-dir is required.", call. = FALSE)
}

study_dir <- normalizePath(args$study_dir, mustWork = TRUE)
cat("Study directory:", study_dir, "\n")
cat("Maximum asset bytes:", format(args$max_asset_bytes, scientific = FALSE), "\n")
cat("Target shard bytes:", format(args$target_shard_bytes, scientific = FALSE), "\n")

publication <- mutt:::build_functional_bundle(
  study_dir = study_dir,
  max_asset_bytes = args$max_asset_bytes,
  target_shard_bytes = args$target_shard_bytes
)

cat("Functional bundle:", publication, "\n")
