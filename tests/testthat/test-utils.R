context("Utility functions")

## Checking for object types---------------------- ------------------------
test_that("Type check reveals wrong admixr type", {
    x <- data.frame()
    class(x) <- "random"
    expect_error(check_type(x, "admixr_result"))
})

test_that("Type check accepts the correct admixr type", {
    x <- data.frame()
    class(x) <- c(class(x), "admixr_result")
    expect_silent(check_type(x, "admixr_result"))
})

# Checking for presence of a sample in an ind file ------------------------

test_that("Fails to find a missing sample", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(admixtools_path(), "data", "allmap"))

  expect_error(check_presence("blah", data))
  expect_error(check_presence(c("blah1", "blah2"), data))
  expect_error(check_presence(c("blah1", "blah2", "French"), data))
})

test_that("Succeeds in finding a sample in an ind file", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(admixtools_path(), "data", "allmap"))

  expect_silent(check_presence("French", data))
  expect_silent(check_presence(c("French", "Neander"), data))
})

# Merging population labels -----------------------------------------------

test_that("Merging of population labels", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(admixtools_path(), "data", "allmap"))

  # list of labels to merge
  merge_list <- list(
    Europe = c("French", "Sardinian", "Czech"),
    WestAfrica = c("Yoruba", "Mende")
  )

  # generate merged ind file using shell utilities
  shell_ind <- tempfile()
  system(stringr::str_c(
    c("sed '",
      as.character(glue::glue("s/{merge_list$Europe}$/Europe/g;")),
      as.character(glue::glue("s/{merge_list$WestAfrica}$/WestAfrica/g;")),
      glue::glue("' {data$ind} > {shell_ind}")),
    collapse = " "
  ))

  # merge labels using
  relabeled <- relabel(data = data,
                       Europe = merge_list$Europe,
                       WestAfrica = merge_list$WestAfrica,
                       outfile = tempfile())

  # compare both generated ind tables
  expect_equivalent(
    readr::read_table(shell_ind, col_names = c("id", "sex", "label")),
    read_ind(relabeled)
  )
})

# SNP counting ------------------------------------------------------------

test_that("SNP counts correspond to numbers from CLI utilities", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(prefix = file.path(admixtools_path(), "convertf", "example"),
                     geno = file.path(admixtools_path(), "convertf", "example.eigenstratgeno"))

  sample_names <- read_ind(data)$id
  shell_counts <- seq_along(sample_names) %>%
    lapply(function(i) {
      paste0("cut -c", i, " ", data$geno, " | sort | uniq -c") %>%
        system(intern = TRUE) %>%
        paste0(collapse = "\n") %>%
        readr::read_table(col_names = c("count", "gt")) %>%
        dplyr::mutate(state = dplyr::case_when(gt == 9 ~ "missing", gt != 9 ~ "present") %>%
                              factor(c("present", "missing"))) %>%
        dplyr::group_by(state) %>%
        dplyr::summarise(sites = sum(count)) %>% tidyr::complete(state, fill = list(sites = 0)) %>%
        dplyr::mutate(id = sample_names[i]) %>%
        tidyr::spread(state, sites)
    }) %>% dplyr::bind_rows()

  expect_equal(count_snps(data)$present, shell_counts$present)
  expect_equal(count_snps(data, missing = TRUE)$missing, shell_counts$missing)
})

# EIGENSTRAT merging--------------------------------------------------------

test_that("Merging produces correct results", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(prefix = file.path(admixtools_path(), "convertf", "example"),
                     geno = file.path(admixtools_path(), "convertf", "example.eigenstratgeno"))

  # read a testing EIGENSTRAT dataset and split it in two parts
  geno <- read_geno(data)
  snp <- read_snp(data)
  ind <- read_ind(data)
  part1 <- list(geno = geno[, 1:2], snp = snp, ind = ind[1:2, ])
  part2 <- list(geno = geno[, 3:ncol(geno)], snp = snp, ind = ind[3:nrow(ind), ])

  # save both parts
  prefix1 <- file.path(tempdir(), "prefix1")
  write_geno(part1$geno, paste0(prefix1, ".geno"))
  write_snp(part1$snp, paste0(prefix1, ".snp"))
  write_ind(part1$ind, paste0(prefix1, ".ind"))
  data1 <- eigenstrat(prefix1)
  prefix2 <- file.path(tempdir(), "prefix2")
  write_geno(part2$geno, paste0(prefix2, ".geno"))
  write_snp(part2$snp, paste0(prefix2, ".snp"))
  write_ind(part2$ind, paste0(prefix2, ".ind"))
  data2 <- eigenstrat(prefix2)

  # merge both parts
  merged_prefix <- file.path(tempdir(), "merged")
  merged <- merge_eigenstrat(merged_prefix, data1, data2)

  # load the merged EIGENSTRAT and compare to the original
  merged_geno <- read_geno(merged)
  merged_snp <- read_snp(merged)
  merged_ind <- read_ind(merged)
  expect_true(all(c(geno == merged_geno, snp == merged_snp, ind == merged_ind)))
})
