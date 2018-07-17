context("Testing the presence of ADMIXTOOLS")

admixtools_present <- function() {
  system("which qpDstat", ignore.stdout = TRUE) == 0
}

path <- admixtools_path()
data_dir <- file.path(path, "data")

test_that("ADMIXTOOLS is present", { expect_true(admixtools_present()) })
test_that("ADMIXTOOLS data is present", { expect_true(dir.exists(data_dir)) })
