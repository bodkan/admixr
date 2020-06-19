context("Examining log outputs")

if (admixtools_present()) {
  path <- admixtools_path()
  data_dir <- file.path(path, "data")
  examples_dir <- file.path(path, "examples")

  setwd(examples_dir)
}

test_that("printlog() generates the same log file as ADMIXTOOLS", {
  skip_on_cran()
  skip_on_os("windows")
  
  data <- eigenstrat(file.path(data_dir, "qpdata"))
  left <- scan(file.path(examples_dir, "left1"), what = "character", quiet = TRUE)
  right <- scan(file.path(examples_dir, "right1"), what = "character", quiet = TRUE) %>%
    stringr::str_subset("^[^#]")

  outdir <- file.path(tempdir(), paste0("log_test", runif(1, 0, 1000000000)))
  dir.create(outdir, showWarnings = FALSE)
  qpAdm_res <- qpAdm(target = left[1], sources = left[-1], outgroups = right, data = data, outdir = outdir)

  loginfo(qpAdm_res, save = TRUE, prefix = "testthat_qpAdm_log", dir = outdir)

  # read the two log files and compare their lines
  log1 <- readLines(list.files(outdir, pattern = "qpAdm__\\d+.log", full.names = TRUE))
  log2 <- readLines(file.path(outdir, "testthat_qpAdm_log.txt"))
  expect_equal(log1, log2)
})



test_that("Requesting log for a missing sample fails", {
  skip_on_cran()
  skip_on_os("windows")
  
  data <- eigenstrat(file.path(data_dir, "qpdata"))
  left <- scan(file.path(examples_dir, "left1"), what = "character", quiet = TRUE)
  right <- scan(file.path(examples_dir, "right1"), what = "character", quiet = TRUE) %>%
    stringr::str_subset("^[^#]")
  
  outdir <- file.path(tempdir(), paste0("log_test", runif(1, 0, 1000000000)))
  dir.create(outdir, showWarnings = FALSE)
  qpAdm_res <- qpAdm(target = left[1], sources = left[-1], outgroups = right, data = data, outdir = outdir)
  
  expect_error(suppressMessages(loginfo(qpAdm_res, target = "missing_sample")))
})



test_that("Specifying target for a non-qpAdm log fails", {
  skip_on_cran()
  skip_on_os("windows")
  
  data <- eigenstrat(file.path(data_dir, "qpdata"))
  left <- scan(file.path(examples_dir, "left1"), what = "character", quiet = TRUE)
  right <- scan(file.path(examples_dir, "right1"), what = "character", quiet = TRUE) %>%
    stringr::str_subset("^[^#]")
  
  outdir <- file.path(tempdir(), paste0("log_test", runif(1, 0, 1000000000)))
  dir.create(outdir, showWarnings = FALSE)
  res <- d(left[1], left[2], left[3], left[4], data = data)
  
  expect_error(suppressMessages(loginfo(res, target = "missing_sample")))
})



test_that("Detect broken log file and inform the user", {
  skip_on_cran()
  skip_on_os("windows")
  
  data <- eigenstrat(file.path(data_dir, "qpdata"))
  truncated_log <- tempfile()
  system(glue::glue('head -n10 {file.path(admixtools_path(), "examples", "qpDstat1.log")} > {truncated_log}'))
  expect_error(capture.output(read_output(truncated_log)))
})
