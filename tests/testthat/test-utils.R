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

  expect_equal(snps_present(prefix)$nsnps, shell_counts$present)
  expect_equal(snps_missing(prefix)$nsnps, shell_counts$missing)
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

test_that("Overlap with no overlapping regions returns an error", {
  prefix <- file.path(admixtools_path(), "convertf", "example")
  bed_file <- tempfile()
  # create a non-sensical BED file that cannot leave any overlapping sites
  data.frame(chrom = c(-1, -1, -1), start = c(1, 10, 20), end = c(2, 11, 21)) %>%
    readr::write_tsv(bed_file, col_names = FALSE)
  # verify that the function fails
  expect_error(subset_sites(prefix = prefix, subset_prefix = "blah", bed_file = bed_file))
})

test_that("Overlap with the same set of sites returns everything", {
  prefix <- file.path(admixtools_path(), "convertf", "example")
  # create a BED file that has the same positions as the original EIGENSTRAT
  bed_file <- tempfile()
  read_snp(paste0(prefix, ".snp")) %>%
    dplyr::mutate(start = pos - 1, end = pos) %>%
    dplyr::select(chrom, start, end) %>%
    readr::write_tsv(bed_file, col_names = FALSE)
  # generate a "subset" based on thad BED file
  subset_prefix <- tempfile()
  subset_sites(prefix = prefix, subset_prefix = subset_prefix, bed_file = bed_file)

  # verify that both EIGENSTRAT datasets are the same
  orig_data <- read_eigenstrat(prefix)
  new_data <- read_eigenstrat(subset_prefix)
  sapply(c("ind", "snp", "geno"), function(i) all(orig_data[[i]] == new_data[[i]])) %>%
    all %>%
    expect_true
})

test_that("Overlap returns a correct number of sites", {
  prefix <- file.path(admixtools_path(), "convertf", "example")
  # resample a BED file a number of times, verifying that we get the correct
  # number of sites after the overlap operation
  orig_data <- read_eigenstrat(prefix)
  successes <- sapply(seq_len(nrow(orig_data$snp)), function(n) {
    # create a BED file that has a subset of positions as the original EIGENSTRAT
    bed_file <- tempfile()
    read_snp(paste0(prefix, ".snp")) %>%
      dplyr::mutate(start = pos - 1, end = pos) %>%
      dplyr::select(chrom, start, end) %>%
      dplyr::sample_n(n) %>%
      dplyr::arrange(chrom, start, end) %>%
      readr::write_tsv(bed_file, col_names = FALSE)
    # generate a "subset" based on thad BED file
    subset_prefix <- tempfile()
    subset_sites(prefix = prefix, subset_prefix = subset_prefix, bed_file = bed_file)

    # verify that both EIGENSTRAT datasets are the same
    new_data <- read_eigenstrat(subset_prefix)
    nrow(new_data$snp) == n
  })
  expect_true(all(successes))
})
