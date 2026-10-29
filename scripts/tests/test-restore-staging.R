test_that('handoff restoration verifies existing assets without rebuilding them', {
  f <- fixture_release('tiny');on.exit(unlink(f$root,recursive=TRUE))
  project <- normalizePath(test_path('..','..'))
  e <- new.env(parent=globalenv())
  withr::with_dir(project,{
    sys.source('R/validation.R',envir=e)
    sys.source('scripts/restore-staging.R',envir=e)
  })
  row <- f$rows$tiny
  jsonlite::write_json(list(schema_version=1L,release_tag='data-v1.0.0',studies=list(row)),
                      file.path(f$root,'assets','studies.candidate.json'),auto_unbox=TRUE)
  staging <- file.path(f$root,'restored')
  expect_output(e$restore_staging(file.path(f$root,'assets'),staging),'verified and restored')
  expect_identical(readRDS(file.path(staging,'tiny','study.rds')),f$objects$tiny)
  expect_error(e$restore_staging(file.path(f$root,'assets'),staging),'overwrite')
  writeLines('corrupt',f$assets$tiny)
  expect_error(e$restore_staging(file.path(f$root,'assets'),file.path(f$root,'failed')),'SHA-256')
  expect_false(dir.exists(file.path(f$root,'failed')))
})
