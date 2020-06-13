context("Examining log outputs")

path <- admixtools_path()
data_dir <- file.path(path, "data")
examples_dir <- file.path(path, "examples")

setwd(examples_dir)

test_that("printlog() generates the same log file as ADMIXTOOLS", {
  skip_on_cran()
  
  data <- eigenstrat(file.path(data_dir, "qpdata"))
  left <- scan(file.path(examples_dir, "left1"), what = "character", quiet = TRUE)
  right <- scan(file.path(examples_dir, "right1"), what = "character", quiet = TRUE) %>%
    stringr::str_subset("^[^#]")

  outdir <- file.path(tempdir(), paste0("log_test", runif(1, 0, 1000000000)))
  dir.create(outdir, showWarnings = FALSE)
  qpAdm_res <- qpAdm(target = left[1], sources = left[-1], outgroups = right, data = data, outdir = outdir)

  printlog(qpAdm_res, save = TRUE, prefix = "testthat_qpAdm_log", dir = outdir)

  # read the two log files and compare their lines
  log1 <- readLines(list.files(outdir, pattern = "qpAdm__\\d+.log", full.names = TRUE))
  log2 <- readLines(file.path(outdir, "testthat_qpAdm_log.txt"))
  expect_equal(log1, log2)
})
