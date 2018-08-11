context("Presence of ADMIXTOOLS on the system")

admixtools_present <- function() {
  system("which qpDstat", ignore.stdout = TRUE) == 0
}

path <- admixtools_path()
data_dir <- file.path(path, "data")

test_that(paste0("ADMIXTOOLS is present\n\n\n\n\n", Sys.getenv("PATH"), "\n\n\n\n",
          system("pwd"), "\n\n\n\n", system("ls")), { expect_true(admixtools_present()) })
test_that("ADMIXTOOLS data is present", { expect_true(dir.exists(data_dir)) })
