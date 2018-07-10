context("Basic ADMIXTOOLS checks")

admixtools_present <- function() {
  system("which qpDstat", ignore.stdout = TRUE) == 0
}

test_that("ADMIXTOOLS is present", { expect_true(admixtools_present()) })
test_that("ADMIXTOOLS data is present", { expect_true(dir.exists(data_path())) })