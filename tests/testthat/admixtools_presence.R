context("Presence of ADMIXTOOLS on the system")

path <- admixtools_path()
data_dir <- file.path(path, "data")

test_that(paste0("ADMIXTOOLS is present", "\n\n\n",
  skip_on_cran()

                 Sys.getenv("PATH"), "\n\n\n",
                 system("ls", intern = TRUE), "\n\n\n",
                 system("pwd", intern = TRUE), "\n\n\n"
                 ), { expect_true(admixtools_present()) })
test_that("ADMIXTOOLS data is present", { expect_true(dir.exists(data_dir)) })
  skip_on_cran()

