context("Basic checks")

admixtools_present <- function() {
  system("which qpDstat", ignore.stdout = TRUE) == 0
}

test_that("ADMIXTOOLS is present", { expect_equal(admixtools_present(), TRUE) })
test_that("ADMIXTOOLS data is present", { expect_equal(dir.exists(data_path()), TRUE) })