#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option("--manifest", type = "character", help = "Candidate studies JSON to verify."),
  make_option("--assets-dir", type = "character", dest = "assets_dir", default = NULL,
              help = "Local release directory. Defaults to the manifest directory."),
  make_option("--remote", action = "store_true", default = FALSE,
              help = "Download each configured URL instead of reading local assets."),
  make_option("--study", type = "character", default = NULL,
              help = "Verify one enabled study instead of every enabled study."),
  make_option("--allow-partial", action = "store_true", dest = "allow_partial", default = FALSE,
              help = "Allow a manifest in which some studies remain disabled."),
  make_option("--promote-to", type = "character", dest = "promote_to", default = NULL,
              help = "After a complete remote check, replace this package manifest.")
)
args <- parse_args(OptionParser(option_list = options))
if (is.null(args$manifest) || !nzchar(args$manifest)) {
  stop("--manifest is required.", call. = FALSE)
}

manifest_path <- normalizePath(args$manifest, mustWork = TRUE)
manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
rows <- as.data.frame(manifest$studies, stringsAsFactors = FALSE)
required <- c("study", "remote_available", "url", "sha256", "size_bytes")
missing_fields <- setdiff(required, names(rows))
if (length(missing_fields)) {
  stop("Manifest is missing fields: ", paste(missing_fields, collapse = ", "), call. = FALSE)
}
if (anyDuplicated(rows$study)) stop("Manifest study IDs must be unique.", call. = FALSE)

enabled <- rows$remote_available %in% TRUE
if (!isTRUE(args$allow_partial) && any(!enabled)) {
  stop("Complete release verification requires every study to be enabled.", call. = FALSE)
}
selected <- if (is.null(args$study)) rows$study[enabled] else args$study
unknown <- setdiff(selected, rows$study[enabled])
if (length(unknown)) {
  stop("Study is absent or not enabled: ", paste(unknown, collapse = ", "), call. = FALSE)
}
if (!length(selected)) stop("No enabled release assets were selected.", call. = FALSE)

assets_dir <- args$assets_dir
if (is.null(assets_dir)) assets_dir <- dirname(manifest_path)
if (!isTRUE(args$remote)) assets_dir <- normalizePath(assets_dir, mustWork = TRUE)

verify_one <- function(study) {
  row <- rows[match(study, rows$study), , drop = FALSE]
  if (!grepl("^[0-9a-fA-F]{64}$", row$sha256[[1L]])) {
    stop("Invalid SHA-256 for `", study, "`.", call. = FALSE)
  }
  if (!is.finite(as.numeric(row$size_bytes[[1L]])) || as.numeric(row$size_bytes[[1L]]) <= 0) {
    stop("Invalid size for `", study, "`.", call. = FALSE)
  }

  temporary <- FALSE
  archive <- if (isTRUE(args$remote)) {
    temporary <- TRUE
    tempfile(paste0(study, "_"), fileext = ".zip")
  } else {
    file.path(assets_dir, paste0(study, ".zip"))
  }
  if (temporary) on.exit(unlink(archive), add = TRUE)
  if (isTRUE(args$remote)) {
    status <- tryCatch(
      utils::download.file(row$url[[1L]], archive, mode = "wb", quiet = FALSE),
      error = identity
    )
    if (inherits(status, "error") || !identical(status, 0L)) {
      detail <- if (inherits(status, "error")) conditionMessage(status) else paste("status", status)
      stop("Remote download failed for `", study, "`: ", detail, call. = FALSE)
    }
  }
  if (!file.exists(archive)) stop("Missing release asset: ", archive, call. = FALSE)
  if (!identical(as.numeric(file.info(archive)$size), as.numeric(row$size_bytes[[1L]]))) {
    stop("Size verification failed for `", study, "`.", call. = FALSE)
  }
  observed <- digest::digest(archive, algo = "sha256", file = TRUE, serialize = FALSE)
  if (!identical(tolower(observed), tolower(row$sha256[[1L]]))) {
    stop("SHA-256 verification failed for `", study, "`.", call. = FALSE)
  }
  members <- utils::unzip(archive, list = TRUE)$Name
  if (!length(members) || any(!startsWith(members, paste0(study, "/")))) {
    stop("Archive layout verification failed for `", study, "`.", call. = FALSE)
  }
  data.frame(study = study, size_bytes = file.info(archive)$size, sha256 = observed)
}

records <- do.call(rbind, lapply(selected, verify_one))
print(records, row.names = FALSE)
cat(sprintf("Verified %d %s release asset(s).\n", nrow(records), if (args$remote) "remote" else "local"))

if (!is.null(args$promote_to)) {
  if (!isTRUE(args$remote)) stop("Manifest promotion requires --remote.", call. = FALSE)
  if (any(!enabled)) stop("Manifest promotion requires a complete release.", call. = FALSE)
  target <- normalizePath(args$promote_to, mustWork = FALSE)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(manifest_path, target, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Could not promote the verified manifest to: ", target, call. = FALSE)
  }
  cat("Promoted verified manifest to ", target, "\n", sep = "")
}
