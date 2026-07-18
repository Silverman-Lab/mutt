make_mutt_fixture <- function(study = "fixture", fail = FALSE) {
  root <- tempfile("mutt_data_")
  parser_root <- tempfile("mutt_parsers_")
  dir.create(file.path(root, study), recursive = TRUE)
  dir.create(file.path(parser_root, study), recursive = TRUE)
  body <- if (fail) {
    sprintf("parse_%s <- function(raw = FALSE, align = FALSE) stop('fixture failure')", study)
  } else {
    sprintf(paste(
      "parse_%s <- function(raw = FALSE, align = FALSE) {",
      "counts <- matrix(c(2, 1), nrow = 1, dimnames = list('S1', c('ASV1', 'ASV2')))",
      "proportions <- counts / rowSums(counts)",
      "scale <- data.frame(load = 3, row.names = 'S1')",
      "metadata <- data.frame(group = 'fixture', row.names = 'S1')",
      "tax <- data.frame(Kingdom = c('Bacteria', 'Bacteria'), row.names = c('ASV1', 'ASV2'))",
      "list(counts = counts, proportions = proportions, scale = scale, metadata = metadata, tax = tax)",
      "}",
      sep = "\n"
    ), study)
  }
  writeLines(body, file.path(parser_root, study, "parse.R"))
  list(root = root, parser_root = parser_root)
}

test_that("mutt is the only exported package function", {
  expect_identical(getNamespaceExports("mutt"), "mutt")
})

test_that("mutt runs a bundled parser without leaking its working directory", {
  fixture <- make_mutt_fixture()
  on.exit(unlink(c(fixture$root, fixture$parser_root), recursive = TRUE), add = TRUE)
  withr::local_options(mutt.parser_root = fixture$parser_root)
  original <- getwd()

  observed <- mutt(studies = "fixture", base_directory = fixture$root)

  expect_s3_class(observed, "mutt_result")
  expect_identical(getwd(), original)
  expect_identical(names(observed), "fixture")
  expect_equal(observed$fixture$counts[1, 1], 2)
  expect_identical(attr(observed, "audit")$status, "success")
})

test_that("parser errors are returned in the structured audit", {
  fixture <- make_mutt_fixture(study = "broken", fail = TRUE)
  on.exit(unlink(c(fixture$root, fixture$parser_root), recursive = TRUE), add = TRUE)
  withr::local_options(mutt.parser_root = fixture$parser_root)
  original <- getwd()

  observed <- NULL
  expect_warning(
    observed <- mutt(studies = "broken", base_directory = fixture$root),
    "1 of 1"
  )

  expect_identical(getwd(), original)
  expect_null(observed$broken)
  expect_identical(attr(observed, "audit")$status, "error")
  expect_match(attr(observed, "audit")$error, "fixture failure")
})

test_that("functional mode and flags are validated before parser execution", {
  expect_error(mutt(rawdata = NA), "rawdata")
  expect_error(mutt(align_samples = 1), "align_samples")
  expect_error(mutt(verbose = NA), "verbose")
  expect_error(mutt(functional = "picrust2"), "FALSE, TRUE")
})
