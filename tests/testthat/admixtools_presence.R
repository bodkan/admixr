context("Presence of ADMIXTOOLS on the system")

test_that("ADMIXTOOLS is present", {
  skip_on_cran()
  expect_true(admixtools_present())
})

test_that("ADMIXTOOLS data is present", {
  skip_on_cran()
  expect_true(dir.exists(file.path(admixtools_path(), "data")))
})
