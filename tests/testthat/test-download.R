test_that('only selected pinned studies download and subsequent requests use cache', {
  f <- fixture_release();on.exit(unlink(f$root,recursive=TRUE))
  withr::local_options(mutt.manifest=f$manifest)
  calls <- character()
  local_mocked_bindings(.mutt_download_file=function(url,destination) {
    id <- sub('[.]zip$','',basename(url));calls <<- c(calls,id)
    file.copy(f$assets[[id]],destination);0L
  },.package='mutt')
  expect_equal(mutt_studies()$study_id,c('tiny','second'))
  one <- mutt('tiny',cache=f$cache)
  expect_identical(one$tiny,f$objects$tiny)
  expect_identical(calls,'tiny')
  expect_identical(mutt('tiny',cache=f$cache),one)
  expect_identical(calls,'tiny')
  both <- mutt(c('tiny','second'),cache=f$cache)
  expect_identical(names(both),c('tiny','second'))
  expect_identical(calls,c('tiny','second'))
  aligned <- mutt('tiny',cache=f$cache,representation='aligned')
  expect_equal(aligned$tiny$metadata$group,c('aligned','aligned'))
  expect_identical(calls,c('tiny','second'))
  mutt_cache_remove('tiny',cache=f$cache)
  expect_true(dir.exists(file.path(f$cache,'data-v1.0.0',paste0('second-',f$rows$second$sha256))))
  mutt('tiny',cache=f$cache)
  expect_identical(calls,c('tiny','second','tiny'))
})

test_that('corrupted and partial downloads never become cache entries', {
  f <- fixture_release('tiny');on.exit(unlink(f$root,recursive=TRUE))
  withr::local_options(mutt.manifest=f$manifest)
  local_mocked_bindings(.mutt_download_file=function(url,destination) {writeLines('incomplete',destination);0L},.package='mutt')
  expect_error(mutt('tiny',cache=f$cache),'SHA-256')
  expect_length(list.files(f$cache,recursive=TRUE,all.files=TRUE),0L)
})

test_that('failed refresh preserves the previous verified copy', {
  f <- fixture_release('tiny');on.exit(unlink(f$root,recursive=TRUE))
  withr::local_options(mutt.manifest=f$manifest)
  fail <- FALSE
  local_mocked_bindings(.mutt_download_file=function(url,destination) {
    if(fail) stop('network unavailable')
    file.copy(f$assets$tiny,destination);0L
  },.package='mutt')
  first <- mutt('tiny',cache=f$cache);fail <- TRUE
  expect_error(mutt('tiny',cache=f$cache,refresh=TRUE),'network unavailable')
  expect_identical(mutt('tiny',cache=f$cache),first)
})

test_that('cache tampering is detected and refresh can repair it', {
  f <- fixture_release('tiny');on.exit(unlink(f$root,recursive=TRUE))
  withr::local_options(mutt.manifest=f$manifest)
  local_mocked_bindings(.mutt_download_file=function(url,destination) {file.copy(f$assets$tiny,destination);0L},.package='mutt')
  mutt('tiny',cache=f$cache)
  target <- file.path(f$cache,'data-v1.0.0',paste0('tiny-',f$rows$tiny$sha256),'study.rds')
  writeLines('tampered',target)
  expect_error(mutt('tiny',cache=f$cache),'checksum')
  expect_identical(mutt('tiny',cache=f$cache,refresh=TRUE)$tiny,f$objects$tiny)
})

test_that('archive path traversal and unlisted members are rejected', {
  f <- fixture_release('tiny');on.exit(unlink(f$root,recursive=TRUE))
  extra <- file.path(f$root,'bad.zip')
  old <- setwd(f$root)
  writeLines('bad','unlisted.txt')
  utils::zip(extra,c('tiny/study.rds','unlisted.txt'),flags='-qX',zip=unname(Sys.which('zip')))
  setwd(old)
  expect_error(mutt:::.mutt_extract_asset(extra,file.path(f$root,'extracted')),'Invalid archive')
})

test_that('manifest rejects mutable, duplicate or incomplete entries', {
  f <- fixture_release('tiny');on.exit(unlink(f$root,recursive=TRUE))
  withr::local_options(mutt.manifest=f$manifest)
  expect_error(mutt(character()),'unique study IDs')
  expect_error(mutt(c('tiny','tiny')),'unique study IDs')
  expect_error(mutt('unknown'),'Unknown')
  expect_error(mutt('tiny',refresh=NA),'refresh')
  manifest <- jsonlite::read_json(f$manifest)
  manifest$studies[[1]]$release_tag <- 'latest'
  jsonlite::write_json(manifest,f$manifest,auto_unbox=TRUE)
  expect_error(mutt_studies(),'Invalid pinned')
})

test_that('draft assets are listed but never advertised as public downloads', {
  f <- fixture_release('tiny');on.exit(unlink(f$root,recursive=TRUE))
  manifest <- jsonlite::read_json(f$manifest);manifest$studies[[1]]$status <- 'staged'
  jsonlite::write_json(manifest,f$manifest,auto_unbox=TRUE)
  withr::local_options(mutt.manifest=f$manifest)
  expect_identical(mutt_studies()$status,'staged')
  expect_error(mutt('tiny',cache=f$cache),'no verified public data release')
})
