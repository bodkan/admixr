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

# SNP counting ------------------------------------------------------------

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

  expect_equal(count_snps(prefix)$present, shell_counts$present)
  expect_equal(count_snps(prefix, missing = TRUE)$missing, shell_counts$missing)
})

# EIGENSTRAT merging--------------------------------------------------------

test_that("Merging produces correct results", {
  prefix <- file.path(admixtools_path(), "convertf", "example")

  # read a testing EIGENSTRAT dataset and split it in two parts
  whole <- read_eigenstrat(prefix)
  part1 <- list(geno = whole$geno[, 1:2], snp = whole$snp, ind = whole$ind[1:2, ])
  part2 <- list(geno = whole$geno[, 3:ncol(whole$geno)], snp = whole$snp, ind = whole$ind[3:nrow(whole$ind), ])

  # save both parts
  prefix1 <- file.path(tempdir(), "prefix1")
  prefix2 <- file.path(tempdir(), "prefix2")
  write_eigenstrat(prefix = prefix1, ind = part1$ind, snp = part1$snp, geno = part1$geno)
  write_eigenstrat(prefix = prefix2, ind = part2$ind, snp = part2$snp, geno = part2$geno)

  # merge both parts
  merged_prefix <- file.path(tempdir(), "merged")
  merge_eigenstrat(prefix = merged_prefix, input1 = prefix1, input2 = prefix2)

  # load the merged EIGENSTRAT and compare to the original
  merged_whole <- read_eigenstrat(merged_prefix)
  expect_true(all(sapply(c("ind", "snp", "geno"), function(i) all(whole[[i]] == merged_whole[[i]]))))
})

# EIGENSTRAT subsetting ---------------------------------------------------

write_bed <- function(snp, bed) {
  read_snp(snp) %>%
    dplyr::mutate(start = pos - 1, end = pos) %>%
    dplyr::select(chrom, start, end) %>%
    readr::write_tsv(bed, col_names = FALSE)
}

test_that("filter_sites correctly handles complete overlap", {
  snp <- paste0(file.path(admixtools_path(), "convertf", "example"), ".snp")
  # create a BED file that has the same positions as the original EIGENSTRAT
  bed <- tempfile()
  write_bed(snp, bed)
  # generate a "subset" based on thad BED file
  output <- tempfile()
  filter_sites(snp, bed, output, include = TRUE)
  # verify that both EIGENSTRAT datasets are the same
  orig_snp <- read_snp(snp)
  output_snp <- read_snp(output)
  expect_equal(orig_snp, output_snp)
})

test_that("filter_sites correctly fails at no overlap", {
  snp <- paste0(file.path(admixtools_path(), "convertf", "example"), ".snp")
  # create a BED file that has the same positions as the original EIGENSTRAT
  bed <- tempfile()
  write_bed(snp, bed)
  # verify that no overlaps leads to error
  expect_error(filter_sites(snp, bed, "blah", include = FALSE))
})

test_that("Overlap returns a correct number of sites", {
  snp_path <- paste0(file.path(admixtools_path(), "convertf", "example"), ".snp")
  snp <- read_snp(snp_path)

  # resample a BED file a number of times, verifying that we get the correct
  # number of sites after the overlap operation
  successes <- sapply(seq_len(nrow(snp)), function(n) {
    # create a BED file that has a subset of sites from the original snp file
    bed <- tempfile()
    snp %>%
      dplyr::mutate(start = pos - 1, end = pos) %>%
      dplyr::select(chrom, start, end) %>%
      dplyr::sample_n(n) %>%
      dplyr::arrange() %>%
      readr::write_tsv(bed, col_names = FALSE)
    # generate a "subset" based on thad BED file
    output <- tempfile()
    filter_sites(snp_path, bed, output, include = TRUE)
    output_snp <- read_snp(output)

    nrow(output_snp) == n
  })
  expect_true(all(successes))
})
