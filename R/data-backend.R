.mutt_data_release <- "2026.07.1"

.mutt_cache_root <- function() {
  override <- getOption("mutt.cache_root", "")
  if (nzchar(override)) return(normalizePath(override, mustWork = FALSE))
  file.path(tools::R_user_dir("mutt", which = "data"), .mutt_data_release, "studies")
}

.mutt_normalize_data_root <- function(path, must_work = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("A data directory must be one nonempty path.", call. = FALSE)
  }
  normalized <- normalizePath(path, mustWork = must_work)
  nested <- file.path(normalized, "studies")
  if (dir.exists(nested)) normalizePath(nested, mustWork = TRUE) else normalized
}

.mutt_local_roots <- function(base_directory = NULL) {
  candidates <- character()
  if (!is.null(base_directory)) {
    candidates <- .mutt_normalize_data_root(base_directory)
  } else {
    configured <- Sys.getenv("MUTT_DATA_DIR", unset = "")
    if (nzchar(configured)) candidates <- c(candidates, .mutt_normalize_data_root(configured))
    candidates <- c(candidates, file.path(getwd(), "studies"), getwd(), .mutt_cache_root())
  }
  candidates <- unique(candidates[dir.exists(candidates)])
  vapply(candidates, normalizePath, character(1), mustWork = TRUE, USE.NAMES = FALSE)
}

.mutt_manifest_path <- function() {
  override <- getOption("mutt.study_manifest", "")
  if (nzchar(override)) return(override)
  installed <- system.file("extdata", "studies.json", package = "mutt")
  if (nzchar(installed)) return(installed)
  local <- file.path(getwd(), "inst", "extdata", "studies.json")
  if (file.exists(local)) return(local)
  ""
}

.mutt_study_manifest <- function() {
  path <- .mutt_manifest_path()
  if (!nzchar(path) || !file.exists(path)) {
    return(data.frame(
      study = .mutt_parser_ids(), remote_available = FALSE,
      url = "", sha256 = "", size_bytes = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  manifest <- jsonlite::read_json(path, simplifyVector = TRUE)
  studies <- as.data.frame(manifest$studies, stringsAsFactors = FALSE)
  required <- c("study", "remote_available", "url", "sha256", "size_bytes")
  missing <- setdiff(required, names(studies))
  if (length(missing)) {
    stop("The MUTT study manifest is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(studies$study)) {
    stop("The MUTT study manifest contains duplicate study IDs.", call. = FALSE)
  }
  enabled <- which(studies$remote_available %in% TRUE)
  if (length(enabled)) {
    valid_sha <- grepl("^[0-9a-fA-F]{64}$", studies$sha256[enabled])
    valid_size <- is.finite(as.numeric(studies$size_bytes[enabled])) &
      as.numeric(studies$size_bytes[enabled]) > 0
    valid_url <- !is.na(studies$url[enabled]) & nzchar(studies$url[enabled])
    if (any(!valid_sha | !valid_size | !valid_url)) {
      stop(
        "Every remotely available MUTT study requires a URL, positive size, and SHA-256.",
        call. = FALSE
      )
    }
  }
  studies
}

.mutt_confirm_all_download <- function(studies, manifest) {
  rows <- match(studies, manifest$study)
  sizes <- suppressWarnings(as.numeric(manifest$size_bytes[rows]))
  total <- sum(sizes, na.rm = TRUE)
  label <- if (is.finite(total) && total > 0) {
    format(structure(total, class = "object_size"), units = "auto")
  } else {
    "an unknown amount of data"
  }
  if (!interactive()) {
    stop(
      "Downloading every MUTT study requires explicit study IDs in noninteractive use.",
      call. = FALSE
    )
  }
  answer <- readline(sprintf(
    "Download %d studies (%s) into the MUTT user cache? [y/N] ",
    length(studies), label
  ))
  if (!tolower(trimws(answer)) %in% c("y", "yes")) {
    stop("Study download cancelled; no files were changed.", call. = FALSE)
  }
  invisible(TRUE)
}

.mutt_sha256 <- function(path) {
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

.mutt_download_study <- function(study, manifest, cache_root) {
  row <- manifest[match(study, manifest$study), , drop = FALSE]
  if (!nrow(row) || !isTRUE(row$remote_available[[1L]]) ||
      !nzchar(row$url[[1L]]) || !nzchar(row$sha256[[1L]])) {
    stop(
      "Study `", study, "` is not available for verified remote download. ",
      "Provide a local Git LFS checkout with `base_directory` or MUTT_DATA_DIR.",
      call. = FALSE
    )
  }

  target <- file.path(cache_root, study)
  if (dir.exists(target)) return(target)
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  lock <- paste0(target, ".lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    stop("Another process is preparing study `", study, "`.", call. = FALSE)
  }
  on.exit(unlink(lock, recursive = TRUE, force = TRUE), add = TRUE)

  staging <- tempfile(paste0(study, "_"), tmpdir = cache_root)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  archive <- file.path(staging, paste0(study, ".zip"))
  status <- tryCatch(
    utils::download.file(row$url[[1L]], archive, mode = "wb", quiet = TRUE),
    error = identity
  )
  if (inherits(status, "error") || !identical(status, 0L)) {
    detail <- if (inherits(status, "error")) conditionMessage(status) else paste("status", status)
    stop("Failed to download study `", study, "`: ", detail, call. = FALSE)
  }
  observed <- .mutt_sha256(archive)
  if (!identical(tolower(observed), tolower(row$sha256[[1L]]))) {
    stop("SHA-256 verification failed for study `", study, "`.", call. = FALSE)
  }

  extracted <- file.path(staging, "extracted")
  dir.create(extracted)
  utils::unzip(archive, exdir = extracted)
  source <- file.path(extracted, study)
  if (!dir.exists(source)) source <- extracted
  if (!file.rename(source, target)) {
    stop("Could not publish the verified study cache for `", study, "`.", call. = FALSE)
  }
  target
}

.mutt_resolve_study_locations <- function(studies, base_directory = NULL,
                                          all_requested = FALSE) {
  roots <- .mutt_local_roots(base_directory)
  locations <- setNames(vector("list", length(studies)), studies)
  for (study in studies) {
    matches <- roots[dir.exists(file.path(roots, study))]
    if (length(matches)) {
      locations[[study]] <- list(
        data_root = matches[[1L]],
        study_dir = file.path(matches[[1L]], study),
        source = if (identical(matches[[1L]], normalizePath(.mutt_cache_root(), mustWork = FALSE))) "cache" else "local"
      )
    }
  }

  missing <- names(locations)[vapply(locations, is.null, logical(1))]
  if (length(missing)) {
    manifest <- .mutt_study_manifest()
    if (all_requested && identical(length(missing), length(studies))) {
      .mutt_confirm_all_download(missing, manifest)
    }
    cache_root <- .mutt_cache_root()
    for (study in missing) {
      .mutt_download_study(study, manifest, cache_root)
      locations[[study]] <- list(
        data_root = cache_root,
        study_dir = file.path(cache_root, study),
        source = "cache"
      )
    }
  }
  locations
}
