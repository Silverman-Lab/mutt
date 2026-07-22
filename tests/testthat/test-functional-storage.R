write_contribution_fixture <- function(path, rows = 30000L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  sample <- sprintf("S%03d", (seq_len(rows) - 1L) %% 12L + 1L)
  function_id <- sprintf("F%04d", (seq_len(rows) - 1L) %% 80L + 1L)
  taxon <- ifelse(seq_len(rows) %% 2L, "P1", "P2")
  value <- sprintf("%.8f", seq_len(rows) / 37)
  table <- data.frame(
    sample = sample,
    `function` = function_id,
    taxon = taxon,
    contribution = value,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.table(
    table, con, sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
  )
  invisible(table)
}

build_functional_publication_fixture <- function(study) {
  branch <- file.path(study, "functional", "picrust2", "asv")
  contribution_paths <- c(
    ec = file.path(branch, "raw", "EC_metagenome_out", "pred_metagenome_contrib.tsv.gz"),
    ko = file.path(branch, "raw", "KO_metagenome_out", "pred_metagenome_contrib.tsv.gz"),
    metacyc_abundance = file.path(branch, "raw", "pathways_out", "path_abun_contrib.tsv.gz")
  )
  tables <- lapply(contribution_paths, write_contribution_fixture)
  dir.create(file.path(branch, "raw", "intermediate"), recursive = TRUE)
  writeLines("not published", file.path(branch, "raw", "intermediate", "large.tmp"))
  dir.create(file.path(branch, "input"), recursive = TRUE)
  writeLines("already stored with study data", file.path(branch, "input", "input.fasta"))
  dir.create(file.path(branch, "logs"), recursive = TRUE)
  writeLines("PICRUSt2 log", file.path(branch, "logs", "stdout.log"))

  descriptors <- Map(
    function(path, type) list(
      relative_path = file.path(
        "raw",
        if (type == "ec") "EC_metagenome_out" else if (type == "ko") "KO_metagenome_out" else "pathways_out",
        basename(path)
      ),
      columns = names(tables[[type]]),
      size_bytes = unname(file.info(path)$size)
    ),
    contribution_paths,
    names(contribution_paths)
  )
  result <- list(
    ec = matrix(1:4, nrow = 2, dimnames = list(c("S001", "S002"), c("E1", "E2"))),
    ko = matrix(5:8, nrow = 2, dimnames = list(c("S001", "S002"), c("K1", "K2"))),
    metacyc_abundance = matrix(1:2, nrow = 2, dimnames = list(c("S001", "S002"), "M1")),
    metacyc_coverage = matrix(1, nrow = 2, dimnames = list(c("S001", "S002"), "M1")),
    stratified = descriptors,
    asv_mapping = data.frame(
      original_feature_id = c("A1", "A2"),
      picrust_id = c("P1", "P2"),
      stringsAsFactors = FALSE
    ),
    taxonomy = data.frame(
      original_feature_id = c("A1", "A2"),
      genus = c("Genus1", "Genus2"),
      stringsAsFactors = FALSE
    ),
    provenance = list(output_directory = branch)
  )
  class(result) <- c("mutt_picrust_branch", "list")
  saveRDS(result, file.path(branch, "result.rds"))
  jsonlite::write_json(
    list(method = "picrust2", branch = "asv", engine_version = "fixture"),
    file.path(branch, "manifest.json"), auto_unbox = TRUE
  )
  jsonlite::write_json(
    data.frame(method = "picrust2", branch = "asv", status = "generated"),
    file.path(study, "functional", "manifest.json"), dataframe = "rows"
  )
  tables
}

test_that("functional publication shards oversized tables and loads transparently", {
  study <- tempfile("mutt_functional_publication_")
  dir.create(study)
  on.exit(unlink(study, recursive = TRUE), add = TRUE)
  original <- build_functional_publication_fixture(study)

  published <- publish_functional_data(
    study,
    max_asset_bytes = 50000,
    target_shard_bytes = 25000
  )
  expect_true(dir.exists(published))
  files <- list.files(published, recursive = TRUE, full.names = TRUE, all.files = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  expect_true(all(unname(file.info(files)$size) <= 50000))
  expect_false(any(grepl("intermediate", files, fixed = TRUE)))
  expect_false(any(grepl("/input/", files, fixed = TRUE)))
  expect_true(file.exists(file.path(published, "publication_manifest.json")))

  loaded <- mutt:::.load_functional_publication(study)
  branch <- loaded$picrust2$asv
  expect_s3_class(branch, "mutt_picrust_branch")
  expect_identical(branch$ec, matrix(1:4, nrow = 2,
                                    dimnames = list(c("S001", "S002"), c("E1", "E2"))))
  expect_identical(branch$stratified$ec$storage, "sharded_tsv_gz")
  expect_error(as.data.frame(branch, type = "ec"), "collect_all")

  selected <- collect_functional(
    branch,
    type = "ec",
    samples = "S001",
    functions = "F0001"
  )
  expect_true(all(selected$sample == "S001"))
  expect_true(all(selected[["function"]] == "F0001"))
  expect_true(all(selected$original_feature_id %in% c("A1", "A2")))
  expect_true("genus" %in% names(selected))

  complete <- collect_functional(branch, type = "ec", collect_all = TRUE)
  expect_equal(nrow(complete), nrow(original$ec))
  expect_identical(sort(unique(complete$sample)), sort(unique(original$ec$sample)))
  expect_identical(complete$sample, original$ec$sample)
  expect_identical(complete[["function"]], original$ec[["function"]])
  expect_identical(complete$picrust_id, original$ec$taxon)
  expect_equal(as.numeric(complete$contribution), as.numeric(original$ec$contribution))

  from_api <- mutt:::.run_study_functional(
    parsed = list(), study_dir = study, mode = "use", tools = list()
  )
  expect_identical(from_api$picrust2$asv$ec, branch$ec)

  first_shard <- file.path(
    branch$provenance$output_directory,
    branch$stratified$ec$shards$relative_path[[1L]]
  )
  unlink(first_shard)
  expect_error(
    collect_functional(branch, type = "ec", collect_all = TRUE),
    "shard is missing"
  )
})

test_that("functional bundle separates core and contribution assets", {
  study <- tempfile("mutt_functional_bundle_")
  original <- build_functional_publication_fixture(study)
  on.exit(unlink(study, recursive = TRUE), add = TRUE)
  withr::local_options(list(
    mutt.functional_core_cache = file.path(study, "core-cache")
  ))

  published <- build_functional_bundle(
    study,
    max_asset_bytes = 50000,
    target_shard_bytes = 25000
  )

  core <- file.path(published, "functional-core.zip")
  manifest_path <- file.path(published, "bundle_manifest.json")
  expect_true(file.exists(core))
  expect_true(file.exists(manifest_path))
  expect_false(file.exists(file.path(published, "publication_index.rds")))
  expect_lte(unname(file.info(core)$size), 50000)

  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  expect_identical(as.integer(manifest$bundle_schema_version), 2L)
  expect_true(is.data.frame(manifest$contributions))
  expect_true(nrow(manifest$contributions) > 3L)
  contribution_files <- file.path(
    published,
    manifest$contributions$publication_path
  )
  expect_true(all(file.exists(contribution_files)))
  expect_true(all(unname(file.info(contribution_files)$size) <= 50000))

  core_members <- utils::unzip(core, list = TRUE)$Name
  expect_true("publication_index.rds" %in% core_members)
  expect_false(any(grepl("stratified", core_members, fixed = TRUE)))
  expect_false(any(grepl("raw/", core_members, fixed = TRUE)))
  expect_false(any(grepl("intermediate", core_members, fixed = TRUE)))

  loaded <- mutt:::.load_functional_publication(study)
  branch <- loaded$picrust2$asv
  expect_s3_class(branch, "mutt_picrust_branch")
  expect_identical(branch$ec, matrix(
    1:4,
    nrow = 2,
    dimnames = list(c("S001", "S002"), c("E1", "E2"))
  ))
  expect_identical(attr(loaded, "publication")$bundle_schema_version, 2L)

  complete <- as.data.frame(branch, type = "ec", collect_all = TRUE)
  expect_equal(nrow(complete), nrow(original$ec))
  expect_identical(complete$sample, original$ec$sample)
  expect_identical(complete[["function"]], original$ec[["function"]])
  expect_identical(complete$picrust_id, original$ec$taxon)
})

test_that("small contribution tables remain single file assets", {
  root <- tempfile("mutt_functional_single_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  source <- file.path(root, "contributions.tsv.gz")
  write_contribution_fixture(source, rows = 20L)
  branch <- file.path(root, "published")
  dir.create(branch)
  published <- mutt:::.functional_publish_tsv(
    source,
    file.path(branch, "ec"),
    branch,
    max_asset_bytes = 100000,
    target_shard_bytes = 50000
  )
  expect_identical(published$descriptor$storage, "single_tsv_gz")
  expect_length(published$descriptor$relative_path, 1L)
  expect_true(file.exists(file.path(branch, published$descriptor$relative_path)))
})

test_that("shard metadata narrows sample, function, and taxon reads", {
  shards <- data.frame(
    relative_path = c("part-1.tsv.gz", "part-2.tsv.gz"),
    stringsAsFactors = FALSE
  )
  shards$sample_ids <- I(list(c("S1", "S2"), c("S3", "S4")))
  shards$function_ids <- I(list(c("F1", "F2"), c("F3", "F4")))
  shards$taxon_ids <- I(list(c("P1", "P2"), c("P3", "P4")))

  expect_identical(
    mutt:::.functional_filter_shards(shards, samples = "S3")$relative_path,
    "part-2.tsv.gz"
  )
  expect_identical(
    mutt:::.functional_filter_shards(shards, functions = "F1")$relative_path,
    "part-1.tsv.gz"
  )
  expect_identical(
    mutt:::.functional_filter_shards(shards, taxa = "P4")$relative_path,
    "part-2.tsv.gz"
  )
  expect_equal(
    nrow(mutt:::.functional_filter_shards(
      shards,
      samples = "S1",
      functions = "F4"
    )),
    0L
  )
})

test_that("oversized result tables are reconstructed from RDS row shards", {
  root <- tempfile("mutt_result_parts_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  set.seed(1)
  result <- list(
    large = matrix(
      stats::runif(12000), nrow = 600,
      dimnames = list(sprintf("S%03d", seq_len(600)), sprintf("F%02d", seq_len(20)))
    ),
    metadata = list(method = "fixture")
  )
  written <- mutt:::.functional_write_result(result, root, max_asset_bytes = 5000)
  expect_identical(written$index$mode, "components")
  restored <- mutt:::.functional_read_result(root, written$index)
  expect_identical(restored, result)
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  expect_true(all(unname(file.info(files)$size) <= 5000))
})
