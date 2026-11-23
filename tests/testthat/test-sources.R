test_that('source metadata exists for every registered study', {
  ids <- mutt_studies()$study_id
  for(id in ids) {
    source <- mutt_source(id)
    expect_identical(source$study_id,id)
    for(item in source$archived_inputs) expect_match(item$sha256,'^[a-f0-9]{64}$')
  }
  expect_error(mutt_source('not-a-study'),'Unknown')
})
