test_that("a verified study archive is downloaded atomically into the cache", {
  study <- "remote_fixture"
  parser_root <- tempfile("mutt_remote_parsers_")
  bundle_root <- tempfile("mutt_remote_bundle_")
  cache_root <- tempfile("mutt_remote_cache_")
  work <- tempfile("mutt_remote_work_")
  dir.create(file.path(parser_root, study), recursive = TRUE)
  dir.create(file.path(bundle_root, study), recursive = TRUE)
  dir.create(work)
  on.exit(unlink(c(parser_root, bundle_root, cache_root, work), recursive = TRUE), add = TRUE)

  writeLines("fixture", file.path(bundle_root, study, "marker.txt"))
  writeLines(sprintf(paste(
    "parse_%s <- function(raw = FALSE, align = FALSE) {",
    "stopifnot(file.exists(file.path('%s', 'marker.txt')))",
    "counts <- matrix(1, 1, 1, dimnames = list('S1', 'ASV1'))",
    "list(counts = counts, proportions = counts,",
    "scale = data.frame(load = 1, row.names = 'S1'),",
    "metadata = data.frame(group = 'fixture', row.names = 'S1'),",
    "tax = data.frame(Kingdom = 'Bacteria', row.names = 'ASV1'))",
    "}", sep = "\n"
  ), study, study), file.path(parser_root, study, "parse.R"))

  archive <- file.path(work, paste0(study, ".zip"))
  withr::with_envvar(
    c(R_ZIPCMD = unname(Sys.which("zip"))),
    withr::with_dir(bundle_root, utils::zip(archive, study, flags = "-rq"))
  )
  manifest <- list(
    schema_version = 1,
    data_release = "fixture",
    studies = list(list(
      study = study,
      parser_version = "fixture",
      remote_available = TRUE,
      url = paste0("file://", archive),
      sha256 = mutt:::.mutt_sha256(archive),
      size_bytes = unname(file.info(archive)$size),
      redistribution_status = "fixture"
    ))
  )
  manifest_path <- file.path(work, "studies.json")
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE)

  withr::local_options(
    mutt.parser_root = parser_root,
    mutt.study_manifest = manifest_path,
    mutt.cache_root = cache_root
  )
  withr::local_envvar(MUTT_DATA_DIR = "")
  withr::local_dir(work)

  observed <- mutt(studies = study)

  expect_identical(attr(observed, "audit")$data_source, "cache")
  expect_true(file.exists(file.path(cache_root, study, "marker.txt")))
  expect_identical(attr(observed, "audit")$status, "success")
})

test_that("remote data with the wrong SHA-256 are rejected", {
  study <- "bad_checksum"
  root <- tempfile("mutt_bad_checksum_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  archive <- file.path(root, "bad.zip")
  writeLines("not a zip", archive)
  manifest <- data.frame(
    study = study,
    remote_available = TRUE,
    url = paste0("file://", archive),
    sha256 = paste(rep("0", 64), collapse = ""),
    size_bytes = file.info(archive)$size,
    stringsAsFactors = FALSE
  )

  expect_error(
    mutt:::.mutt_download_study(study, manifest, file.path(root, "cache")),
    "SHA-256 verification failed"
  )
})

test_that("enabled remote entries require complete integrity metadata", {
  root <- tempfile("mutt_incomplete_manifest_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  manifest_path <- file.path(root, "studies.json")
  jsonlite::write_json(
    list(
      schema_version = 1,
      data_release = "fixture",
      studies = list(list(
        study = "incomplete",
        remote_available = TRUE,
        url = "https://example.org/incomplete.zip",
        sha256 = "",
        size_bytes = 0
      ))
    ),
    manifest_path,
    auto_unbox = TRUE
  )
  withr::local_options(mutt.study_manifest = manifest_path)

  expect_error(
    mutt:::.mutt_study_manifest(),
    "requires a URL, positive size, and SHA-256"
  )
})
