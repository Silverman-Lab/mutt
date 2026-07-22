.functional_bundle_version <- 2L
.functional_core_archive <- "functional-core.zip"
.functional_bundle_manifest <- "bundle_manifest.json"

.functional_full_asset_path <- function(branch, relative_path) {
  gsub(
    "\\\\",
    "/",
    file.path(as.character(branch), as.character(relative_path)),
    fixed = FALSE
  )
}

.functional_validate_archive_members <- function(members) {
  normalized <- gsub("\\\\", "/", members)
  unsafe <- startsWith(normalized, "/") |
    grepl("^[A-Za-z]:/", normalized) |
    vapply(
      strsplit(normalized, "/", fixed = TRUE),
      function(parts) any(parts == ".."),
      logical(1)
    )
  if (any(unsafe)) {
    stop("Functional core archive contains an unsafe path.", call. = FALSE)
  }
  invisible(TRUE)
}

.functional_copy_relative_files <- function(source_root, destination_root, files) {
  for (relative in files) {
    source <- file.path(source_root, relative)
    destination <- file.path(destination_root, relative)
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source, destination, overwrite = TRUE, copy.mode = TRUE)) {
      stop("Could not copy functional bundle asset: ", source, call. = FALSE)
    }
  }
  invisible(TRUE)
}

.functional_pack_publication <- function(study_dir, max_asset_bytes) {
  publication_root <- file.path(study_dir, "functional_data")
  index_path <- file.path(publication_root, "publication_index.rds")
  if (!file.exists(index_path)) {
    stop("Loose functional publication index was not created.", call. = FALSE)
  }

  index <- readRDS(index_path)
  assets <- index$assets
  if (!is.data.frame(assets) ||
      !all(c("branch", "relative_path", "role") %in% names(assets))) {
    stop("Functional publication has an invalid asset index.", call. = FALSE)
  }

  asset_paths <- .functional_full_asset_path(assets$branch, assets$relative_path)
  contribution_rows <- grepl("^stratified_", as.character(assets$role))
  contribution_paths <- unique(asset_paths[contribution_rows])

  all_files <- list.files(
    publication_root,
    recursive = TRUE,
    all.files = TRUE,
    no.. = TRUE,
    include.dirs = FALSE
  )
  all_files <- gsub("\\\\", "/", all_files)
  core_files <- setdiff(all_files, contribution_paths)
  if (!length(core_files)) {
    stop("Functional core contains no files.", call. = FALSE)
  }

  staging <- tempfile("functional-bundle-staging-", tmpdir = study_dir)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)

  .functional_copy_relative_files(publication_root, staging, contribution_paths)

  core_archive <- file.path(staging, .functional_core_archive)
  zip_command <- unname(Sys.which("zip"))
  if (!nzchar(zip_command)) {
    stop("The `zip` command is required to build a functional bundle.", call. = FALSE)
  }
  status <- withr::with_envvar(
    c(R_ZIPCMD = zip_command),
    withr::with_dir(
      publication_root,
      utils::zip(core_archive, core_files, flags = "-qX9")
    )
  )
  if (!identical(status, 0L) || !file.exists(core_archive)) {
    stop("Could not create functional core archive.", call. = FALSE)
  }
  core_size <- unname(file.info(core_archive)$size)
  if (!is.finite(core_size) || core_size > max_asset_bytes) {
    stop(
      "Functional core archive exceeds per-file publication limit: ",
      core_size,
      " bytes.",
      call. = FALSE
    )
  }

  contribution_manifest <- assets[contribution_rows, , drop = FALSE]
  contribution_manifest$publication_path <- asset_paths[contribution_rows]
  contribution_manifest <- contribution_manifest[, c(
    "publication_path", "branch", "relative_path", "role", "source",
    "size_bytes", "sha256", "rows"
  )]

  bundle <- list(
    bundle_schema_version = .functional_bundle_version,
    publication_schema_version = index$schema_version,
    created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    core = list(
      relative_path = .functional_core_archive,
      size_bytes = core_size,
      sha256 = .functional_sha256(core_archive)
    ),
    contributions = contribution_manifest
  )
  jsonlite::write_json(
    bundle,
    file.path(staging, .functional_bundle_manifest),
    dataframe = "rows",
    auto_unbox = TRUE,
    pretty = TRUE,
    na = "null"
  )

  published_files <- list.files(
    staging,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  published_files <- published_files[file.info(published_files)$isdir %in% FALSE]
  oversized <- published_files[unname(file.info(published_files)$size) > max_asset_bytes]
  if (length(oversized)) {
    stop(
      "Functional bundle contains oversized asset(s): ",
      paste(oversized, collapse = ", "),
      call. = FALSE
    )
  }

  backup <- paste0(publication_root, ".loose")
  if (dir.exists(backup)) unlink(backup, recursive = TRUE, force = TRUE)
  if (!file.rename(publication_root, backup)) {
    stop("Could not protect loose functional publication.", call. = FALSE)
  }
  if (!file.rename(staging, publication_root)) {
    file.rename(backup, publication_root)
    stop("Could not publish functional bundle atomically.", call. = FALSE)
  }
  unlink(backup, recursive = TRUE, force = TRUE)
  invisible(publication_root)
}

#' Build a size-safe functional-data bundle
#'
#' Internal maintainer helper used by the HPC finalizer. It first constructs
#' the validated loose publication, then stores metadata and compact results in
#' `functional-core.zip` while retaining EC, KO, and MetaCyc contribution
#' assets separately. Oversized contribution tables remain sharded.
#'
#' @inheritParams publish_functional_data
#' @return Path to the published `functional_data/` bundle, invisibly.
#' @keywords internal
build_functional_bundle <- function(
  study_dir,
  max_asset_bytes = .functional_default_max_asset_bytes,
  target_shard_bytes = min(.functional_default_shard_bytes, max_asset_bytes * 0.75)
) {
  max_asset_bytes <- .functional_validate_bytes(max_asset_bytes, "max_asset_bytes")
  target_shard_bytes <- .functional_validate_bytes(target_shard_bytes, "target_shard_bytes")
  publish_functional_data(
    study_dir = study_dir,
    max_asset_bytes = max_asset_bytes,
    target_shard_bytes = target_shard_bytes
  )
  .functional_pack_publication(study_dir, max_asset_bytes)
}

.functional_read_bundle_manifest <- function(publication_root) {
  path <- file.path(publication_root, .functional_bundle_manifest)
  if (!file.exists(path)) return(NULL)
  manifest <- jsonlite::read_json(path, simplifyVector = TRUE)
  if (!identical(as.integer(manifest$bundle_schema_version), .functional_bundle_version)) {
    stop("Unsupported functional bundle schema version.", call. = FALSE)
  }
  manifest
}

.functional_extract_core <- function(publication_root, bundle) {
  core_archive <- file.path(publication_root, bundle$core$relative_path)
  if (!file.exists(core_archive)) {
    stop("Functional core archive is missing: ", core_archive, call. = FALSE)
  }
  observed_size <- unname(file.info(core_archive)$size)
  if (!identical(as.numeric(observed_size), as.numeric(bundle$core$size_bytes))) {
    stop("Functional core archive size does not match its manifest.", call. = FALSE)
  }

  cache_parent <- getOption(
    "mutt.functional_core_cache",
    file.path(tools::R_user_dir("mutt", "cache"), "functional-core")
  )
  cache_key <- as.character(bundle$core$sha256)
  cache_root <- file.path(cache_parent, cache_key)
  complete <- file.path(cache_root, ".complete")
  if (file.exists(complete)) return(cache_root)

  observed_hash <- .functional_sha256(core_archive)
  if (!identical(tolower(observed_hash), tolower(cache_key))) {
    stop("Functional core archive checksum does not match its manifest.", call. = FALSE)
  }

  dir.create(cache_parent, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cache_parent)) {
    stop("Could not create functional core cache: ", cache_parent, call. = FALSE)
  }
  staging <- tempfile("functional-core-", tmpdir = cache_parent)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  on.exit(if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)

  members <- utils::unzip(core_archive, list = TRUE)$Name
  .functional_validate_archive_members(members)
  utils::unzip(core_archive, exdir = staging)
  if (!file.exists(file.path(staging, "publication_index.rds"))) {
    stop("Functional core archive has no publication index.", call. = FALSE)
  }
  writeLines(cache_key, file.path(staging, ".complete"))

  if (!file.rename(staging, cache_root)) {
    if (!file.exists(complete)) {
      stop("Could not install extracted functional core cache.", call. = FALSE)
    }
  }
  cache_root
}

.functional_publication_roots <- function(study_dir) {
  publication_root <- file.path(study_dir, "functional_data")
  bundle <- .functional_read_bundle_manifest(publication_root)
  if (is.null(bundle)) {
    return(list(core_root = publication_root, asset_root = publication_root, bundle = NULL))
  }
  list(
    core_root = .functional_extract_core(publication_root, bundle),
    asset_root = publication_root,
    bundle = bundle
  )
}
