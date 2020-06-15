# This script simply runs the test suite that is part of the ADMIXTOOLS
# package. It is normally executed via Perl script in examples/mklog.

context("Built-in ADMIXTOOLS test suite")

if (admixtools_present()) {
  path <- admixtools_path()
  data_dir <- file.path(path, "data")
  examples_dir <- file.path(path, "examples")
}

test_that("qpDstat test passes (4 population input version)", {
  skip_on_cran()
  skip_on_os("windows")

  setwd(examples_dir)
  expect_silent(run_cmd("qpDstat",
                        par_file = file.path(examples_dir, "parqpDstat1"),
                        log_file = file.path(examples_dir, "test_qpDstat1.log")))
})

test_that("qp3Pop test passes", {
  skip_on_cran()
  skip_on_os("windows")

  setwd(examples_dir)
  expect_silent(run_cmd("qp3Pop",
                        par_file = file.path(examples_dir, "parqp3Pop"),
                        log_file = file.path(examples_dir, "test_qp3Pop.log")))
})

test_that("qpF4ratio test passes", {
  skip_on_cran()
  skip_on_os("windows")

  setwd(examples_dir)
  expect_silent(run_cmd("qpF4ratio",
                        par_file = file.path(examples_dir, "parqpF4ratio"),
                        log_file = file.path(examples_dir, "test_qpF4ratio.log")))
})

test_that("qpAdm test passes", {
  skip_on_cran()
  skip_on_os("windows")

  setwd(examples_dir)
  expect_silent(run_cmd("qpAdm",
                        par_file = file.path(examples_dir, "parqpAdm"),
                        log_file = file.path(examples_dir, "test_qpAdm.log")))
})

test_that("qpWave test passes", {
  skip_on_cran()
  skip_on_os("windows")

  setwd(examples_dir)
  expect_silent(run_cmd("qpWave",
                        par_file = file.path(examples_dir, "parqpWave"),
                        log_file = file.path(examples_dir, "test_qpWave.log")))
})
