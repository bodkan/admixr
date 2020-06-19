# These tests run the same calculations as those performed by ADMIXTOOLS'
# built in tests (see test-01_admixtools_cmds.R), but they do so using
# admixr's wrapper functions. The results of both tests (i.e. contents of
# output log files) are then compared to each other.

context("Wrapper functionality")

if (admixtools_present()) {
  path <- admixtools_path()
  data_dir <- file.path(path, "data")
  examples_dir <- file.path(path, "examples")

  setwd(examples_dir)
}

# f3 ----------------------------------------------------------------------

test_that("qp3Pop wrapper produces correct results", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(data_dir, "allmap"))
  pops <- read_pops(file.path(examples_dir, "list_qp3Pop"), c("A", "B", "C"))
  expect_equivalent(
    f3(A = pops$A, B = pops$B, C = pops$C, data = data),
    read_output(file.path(examples_dir, "test_qp3Pop.log"))
  )
})

# D -----------------------------------------------------------------------

test_that("qpDstat wrapper produces correct results (4 population input version)", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(data_dir, "allmap"))
  pops <- read_pops(file.path(examples_dir, "list_qpDstat1"), c("W", "X", "Y", "Z"))
  expect_equivalent(
    dplyr::select(d(W = pops$W, X = pops$X, Y = pops$Y, Z = pops$Z, data = data), -stderr),
    read_output(file.path(examples_dir, "test_qpDstat1.log"))
  )
})

# f4-ratio ----------------------------------------------------------------

test_that("qpF4ratio wrapper produces correct results", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(data_dir, "allmap"))
  pops <- readLines(file.path(examples_dir, "list_qpF4ratio")) %>%
    stringr::str_replace_all(" :+ ", " ") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_c(collapse = "\n") %>%
    readr::read_delim(
      delim = " ",
      col_names = c("A", "O", "X", "C", "_", "__", "B", "___")
    ) %>%
    dplyr::select(X, A, B, C, O)
  expect_equivalent(
    dplyr::bind_rows(lapply(seq_len(nrow(pops)), function(i) {
      f4ratio(
        X = pops$X[i], A = pops$A[i], B = pops$B[i], C = pops$C[i], O = pops$O[i],
        data = data
      )
    })),
    read_output(file.path(examples_dir, "test_qpF4ratio.log"))
  )
})

# qpAdm -------------------------------------------------------------------

test_that("qpAdm wrapper produces correct results", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(data_dir, "qpdata"))
  left <- scan(file.path(examples_dir, "left1"), what = "character", quiet = TRUE)
  right <- scan(file.path(examples_dir, "right1"), what = "character", quiet = TRUE) %>%
    stringr::str_subset("^[^#]")
  props1 <- qpAdm(target = left[1], sources = left[-1], outgroups = right, data = data, params=NULL)$proportions
  props2 <- read_output(file.path(examples_dir, "test_qpAdm.log"))$proportions
  expect_equal(props1, props2)
})


test_that("qpAdm with a single source produces NULL subsets dataframe", {
  skip_on_cran()
  skip_on_os("windows")
  
  data <- eigenstrat(file.path(data_dir, "qpdata"))
  left <- scan(file.path(examples_dir, "left1"), what = "character", quiet = TRUE)[1:2]
  right <- scan(file.path(examples_dir, "right1"), what = "character", quiet = TRUE) %>%
    stringr::str_subset("^[^#]")
  result <- qpAdm(target = left[1], sources = left[-1], outgroups = right, data = data)
  expect_true(is.null(result$subsets))
})


# qpWave ------------------------------------------------------------------

test_that("qpWave wrapper produces correct results", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(data_dir, "qpdata"))
  left <- scan(file.path(examples_dir, "left1"), what = "character", quiet = TRUE)
  right <- scan(file.path(examples_dir, "right1"), what = "character", quiet = TRUE) %>%
    stringr::str_subset("^[^#]")
  expect_equivalent(
    qpWave(left = left, right = right, data = data),
    read_output(file.path(examples_dir, "test_qpWave.log"))
  )
})


# test that the output object is of the right class ----------------------
test_that("Output object is of the right class", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(data_dir, "allmap"))
  pops <- read_pops(file.path(examples_dir, "list_qp3Pop"), c("A", "B", "C"))
  res <- f3(A = pops$A, B = pops$B, C = pops$C, data = data)
  expect_true(is(res, "admixr_result"))
})

# test that the output object gets proper annotation ----------------------
test_that("Output object is properly annotated with the command name", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(data_dir, "allmap"))
  pops <- read_pops(file.path(examples_dir, "list_qp3Pop"), c("A", "B", "C"))
  res <- f3(A = pops$A, B = pops$B, C = pops$C, data = data)
  expect_true(attr(res, "command") == "qp3Pop")
})

# test that the output object contains the whole log ----------------------
test_that("Output object carries with it the whole log output", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(data_dir, "allmap"))
  pops <- read_pops(file.path(examples_dir, "list_qp3Pop"), c("A", "B", "C"))
  outdir <- file.path(tempdir(), paste0("log_test", runif(1, 0, 10000000)))
  dir.create(outdir)
  res <- f3(A = pops$A, B = pops$B, C = pops$C, data = data, outdir = outdir)
  expect_true(all(attr(res, "log_output") == readLines(list.files(outdir, pattern = "*.log$", full.names = TRUE))))
})

