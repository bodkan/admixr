context("Utility functions")

# Checking for presence of a sample in an ind file ------------------------

test_that("Fails to find a missing sample", {
  prefix <- file.path(admixtools_path(), "data", "allmap")
  expect_error(check_presence("blah", ind = paste0(prefix, ".ind")))
  expect_error(check_presence("blah", prefix))
  expect_error(check_presence(c("blah1", "blah2"), prefix))
  expect_error(check_presence(c("blah1", "blah2", "French"), prefix))
})

test_that("Succeeds in finding a sample in an ind file", {
  prefix <- file.path(admixtools_path(), "data", "allmap")
  expect_silent(check_presence("French", ind = paste0(prefix, ".ind")))
  expect_silent(check_presence("French", prefix))
  expect_silent(check_presence(c("French", "Neander"), prefix))
})

# Merging population labels -----------------------------------------------

test_that("Merging of population labels", {
  prefix <- file.path(admixtools_path(), "data", "allmap")

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
      glue::glue("' {prefix}.ind > {shell_ind}")),
    collapse = " "
  ))

  # merge labels using
  admixr_ind <- tempfile()
  merge_labels(ind = paste0(prefix, ".ind"),
               modified_ind = admixr_ind,
               labels = merge_list)

  # compare both generated ind tables
  expect_equal(read_ind(shell_ind), read_ind(admixr_ind))
})

# SNP counting -----------------------------------------------

test_that("SNP counts correspond to numbers from CLI utilities", {
  prefix <- file.path(admixtools_path(), "convertf", "example")
  sample_names <- read_ind(paste0(prefix, ".ind"))$id
  shell_counts <- seq_along(sample_names) %>%
    lapply(function(i) {
      paste0("cut -c", i, " ", prefix, ".geno | sort | uniq -c") %>%
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

  expect_equal(snps_present(prefix)$nsnps, shell_counts$present)
  expect_equal(snps_missing(prefix)$nsnps, shell_counts$missing)
})
