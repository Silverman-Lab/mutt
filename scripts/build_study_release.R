#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option("--study", type = "character",
              help = "One study ID to bundle. Mutually exclusive with --all."),
  make_option("--all", action = "store_true", default = FALSE,
              help = "Bundle every registered study as one data release."),
  make_option("--check-only", action = "store_true", dest = "check_only", default = FALSE,
              help = "Validate release readiness without writing archives."),
  make_option("--data-root", type = "character", dest = "data_root", default = "studies",
              help = "Directory containing study folders [default: %default]."),
  make_option("--manifest", type = "character", default = "inst/extdata/studies.json",
              help = "Package study manifest [default: %default]."),
  make_option("--output-dir", type = "character", dest = "output_dir", default = NULL,
              help = "New release staging directory [default: release/data-v<release>].")
)
args <- parse_args(OptionParser(option_list = options))

if (isTRUE(args$all) == !is.null(args$study)) {
  stop("Supply exactly one of --study STUDY_ID or --all.", call. = FALSE)
}

manifest_path <- normalizePath(args$manifest, mustWork = TRUE)
manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
rows <- as.data.frame(manifest$studies, stringsAsFactors = FALSE)
required <- c(
  "study", "parser_file", "remote_available", "url", "sha256",
  "size_bytes", "redistribution_status"
)
missing_fields <- setdiff(required, names(rows))
if (length(missing_fields)) {
  stop("Manifest is missing fields: ", paste(missing_fields, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(rows$study)) stop("Manifest study IDs must be unique.", call. = FALSE)

selected <- if (isTRUE(args$all)) rows$study else args$study
unknown <- setdiff(selected, rows$study)
if (length(unknown)) stop("Unknown study: ", paste(unknown, collapse = ", "), call. = FALSE)

data_root <- normalizePath(args$data_root, mustWork = TRUE)
parser_root <- normalizePath(file.path(dirname(manifest_path), "..", "parsers"), mustWork = TRUE)
selected_rows <- rows[match(selected, rows$study), , drop = FALSE]
readiness <- data.frame(
  study = selected,
  study_directory = dir.exists(file.path(data_root, selected)),
  parser = file.exists(file.path(parser_root, selected_rows$parser_file)),
  redistribution = selected_rows$redistribution_status == "verified",
  stringsAsFactors = FALSE
)
readiness$ready <- readiness$study_directory & readiness$parser & readiness$redistribution

cat(sprintf(
  "Data release %s: %d selected; %d ready; %d blocked.\n",
  manifest$data_release, nrow(readiness), sum(readiness$ready), sum(!readiness$ready)
))
if (any(!readiness$ready)) print(readiness[!readiness$ready, ], row.names = FALSE)

if (isTRUE(args$check_only)) {
  if (any(!readiness$ready)) quit(save = "no", status = 2L)
  quit(save = "no", status = 0L)
}
if (any(!readiness$ready)) {
  stop(
    "Release build stopped before writing archives. Resolve every blocked row or build one verified study.",
    call. = FALSE
  )
}

output_dir <- args$output_dir
if (is.null(output_dir)) {
  output_dir <- file.path("release", paste0("data-v", manifest$data_release))
}
output_dir <- normalizePath(output_dir, mustWork = FALSE)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Refusing to write into a nonempty release directory: ", output_dir, call. = FALSE)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

zip_command <- unname(Sys.which("zip"))
if (!nzchar(zip_command)) stop("The `zip` command is required.", call. = FALSE)

build_one <- function(study) {
  study_dir <- file.path(data_root, study)
  files <- list.files(study_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  excluded <- grepl("(^|/)(functional|orientation_repair[.]staging)(/|$)", files) |
    grepl("(^|/)([.]DS_Store|[.]Rhistory)$", files)
  files <- files[!excluded]
  if (!length(files)) stop("No release files remain for `", study, "`.", call. = FALSE)

  archive <- file.path(output_dir, paste0(study, ".zip"))
  relative <- file.path(study, files)
  status <- withr::with_envvar(
    c(R_ZIPCMD = zip_command),
    withr::with_dir(data_root, utils::zip(archive, relative, flags = "-qX9"))
  )
  if (!identical(status, 0L) || !file.exists(archive)) {
    stop("Archive creation failed for `", study, "`.", call. = FALSE)
  }

  members <- utils::unzip(archive, list = TRUE)$Name
  if (!length(members) || any(!startsWith(members, paste0(study, "/")))) {
    stop("Archive layout validation failed for `", study, "`.", call. = FALSE)
  }
  data.frame(
    study = study,
    archive = normalizePath(archive),
    asset = basename(archive),
    size_bytes = unname(file.info(archive)$size),
    sha256 = digest::digest(archive, algo = "sha256", file = TRUE, serialize = FALSE),
    stringsAsFactors = FALSE
  )
}

records <- do.call(rbind, lapply(selected, build_one))
candidate <- manifest
candidate_rows <- as.data.frame(candidate$studies, stringsAsFactors = FALSE)
for (i in seq_len(nrow(records))) {
  j <- match(records$study[[i]], candidate_rows$study)
  candidate_rows$remote_available[[j]] <- TRUE
  candidate_rows$size_bytes[[j]] <- records$size_bytes[[i]]
  candidate_rows$sha256[[j]] <- records$sha256[[i]]
}
candidate$studies <- lapply(seq_len(nrow(candidate_rows)), function(i) {
  as.list(candidate_rows[i, , drop = FALSE])
})

jsonlite::write_json(
  list(
    data_release = manifest$data_release,
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    assets = lapply(seq_len(nrow(records)), function(i) as.list(records[i, , drop = FALSE]))
  ),
  file.path(output_dir, "release-assets.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
jsonlite::write_json(
  candidate,
  file.path(output_dir, "studies.candidate.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

print(records, row.names = FALSE)
cat("\nBuilt assets are local staging files, not published remote assets.\n")
cat("After upload, verify every URL and checksum before replacing the package manifest.\n")
