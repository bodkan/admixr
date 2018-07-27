# These tests run the same calculations as those performed by ADMIXTOOLS'
# built in tests (see test-01_admixtools_cmds.R), but they do so using
# admixr's wrapper functions. The results of both tests (i.e. contents of
# output log files) are then compared to each other.

context("Testing the wrapper functions")

path <- admixtools_path()
data_dir <- file.path(path, "data")
examples_dir <- file.path(path, "examples")

setwd(examples_dir)

read_pops <- function(filename, columns) {
  df <- setNames(read.table(filename, stringsAsFactors = FALSE), columns)
  lapply(df, function(col) unique(col))
}

test_that("qp3Pop wrapper produces correct results", {
  prefix <- paste0(data_dir, "/allmap")
  pops <- read_pops(file.path(examples_dir, "list_qp3Pop"), c("A", "B", "C"))
  expect_equal(
    qp3Pop(A = pops$A, B = pops$B, C = pops$C, prefix = prefix),
    read_output(file.path(examples_dir, "test_qp3Pop.log"))
  )
})

test_that("qpDstat wrapper produces correct results (4 population input version)", {
  prefix <- paste0(data_dir, "/allmap")
  pops <- read_pops(file.path(examples_dir, "list_qpDstat1"), c("W", "X", "Y", "Z"))
  expect_equal(
    dplyr::select(qpDstat(W = pops$W, X = pops$X, Y = pops$Y, Z = pops$Z, prefix = prefix), -stderr),
    read_output(file.path(examples_dir, "test_qpDstat1.log"))
  )
})

test_that("qpF4ratio wrapper produces correct results", {
  prefix <- paste0(data_dir, "/allmap")
  pops <- readLines(file.path(examples_dir, "list_qpF4ratio")) %>%
    stringr::str_replace_all(" :+ ", " ") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_c(collapse = "\n") %>%
    readr::read_delim(
      delim = " ",
      col_names = c("A", "O", "X", "C", "_", "__", "B", "___")
    ) %>%
    dplyr::select(X, A, B, C, O)
  expect_equal(
    dplyr::bind_rows(lapply(seq_len(nrow(pops)), function(i) {
      qpF4ratio(
        X = pops$X[i], A = pops$A[i], B = pops$B[i], C = pops$C[i], O = pops$O[i],
        prefix = prefix
      )
    })),
    read_output(file.path(examples_dir, "test_qpF4ratio.log"))
  )
})

test_that("qpAdm wrapper produces correct results", {
  prefix <- paste0(data_dir, "/qpdata")
  left <- scan(file.path(examples_dir, "left1"), what = "character", quiet = TRUE)
  right <- scan(file.path(examples_dir, "right1"), what = "character", quiet = TRUE) %>%
    stringr::str_subset("^[^#]")
  expect_equal(
    qpAdm(target = left[1], source = left[-1], outgroup = right, prefix = prefix),
    read_output(file.path(examples_dir, "test_qpAdm.log"))
  )
})
