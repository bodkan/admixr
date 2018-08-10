context("Utility functions")

# Checking for presence of a sample in an ind file ------------------------

test_that("Fails to find a missing sample", {
  prefix <- file.path(admixtools_path(), "data", "allmap")
  expect_error(check_presence("blah", ind = paste0(prefix, ".ind")))
  expect_error(check_presence("blah", prefix))
  expect_error(check_presence(c("blah1", "blah2"), prefix))
  expect_error(check_presence(c("blah1", "blah2", "French"), prefix))
})

test_that("Succeeds in finding a sample in an ind file", {
  prefix <- file.path(admixtools_path(), "data", "allmap")
  expect_silent(check_presence("French", ind = paste0(prefix, ".ind")))
  expect_silent(check_presence("French", prefix))
  expect_silent(check_presence(c("French", "Neander"), prefix))
})
