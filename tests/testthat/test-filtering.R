context("Filtering of sites")

read_snp_file <- function(path) {
  readr::read_table2(
    path,
    col_types = "ccdicc",
    col_names = c("id", "chrom", "gen", "pos", "ref", "alt")
  )
}

snp_to_bed <- function(snp, bed) {
  read_snp_file(snp) %>%
        dplyr::mutate(start = pos - 1, end = pos) %>%
        dplyr::select(chrom, start, end) %>%
        readr::write_tsv(bed, col_names = FALSE)
}

test_that("Potential aDNA SNPs are correctly removed", {
  data <- eigenstrat(file.path(admixtools_path(), "data/allmap"))
  output <- tempfile()
  system(
    sprintf("awk '($5 == \"C\" && $6 == \"T\") ||
                  ($5 == \"T\" && $6 == \"C\") ||
                  ($5 == \"G\" && $6 == \"A\") ||
                   ($5 == \"A\" && $6 == \"G\")' %s > %s", data$snp, output)
  )
  n_shell <- nrow(read_snp_file(output))
  n_admixr <- nrow(read_snp_file(remove_transitions(data)$exclude))
  expect_equal(n_shell, n_admixr)
})


orig_data <- eigenstrat(file.path(admixtools_path(), "convertf", "example"))

# ADMIXTOOLS example data is broken and it's first SNP has a position 0,
# although snp files have to be 0-based - let's remove the first SNP entirely
data <- orig_data
data$snp <- tempfile()
data$geno <- tempfile()
system(sprintf("tail -n+2 %s > %s", orig_data$snp, data$snp))
system(sprintf("tail -n+2 %s > %s", orig_data$geno, data$geno))

test_that("filter_sites correctly fails at no overlap", {
    # create a BED file that has the same positions as the original EIGENSTRAT
    bed <- tempfile()
    snp_to_bed(data$snp, bed)

    # verify that no overlaps leads to error
    expect_error(filter_bed(data, bed, remove = TRUE))
})

test_that("filter_sites correctly handles complete overlap", {
    # create a BED file that has the same positions as the original EIGENSTRAT
    bed <- tempfile()
    snp_to_bed(data$snp, bed)

    # generate a "subset" based on that BED file
    output <- tempfile()
    new_data <- filter_bed(data = data, bed = bed)

    # verify that both EIGENSTRAT datasets are the same
    expect_true(nrow(read_snp_file(new_data$exclude)) == 0)
})

test_that("Overlap returns a correct number of sites", {
    snp <- read_snp_file(data$snp)

    # resample a BED file a number of times, verifying that we get the correct
    # number of sites after the overlap operation
    successes <- sapply(seq_len(nrow(snp) - 1), function(n) {
        # create a BED file that has a subset of sites from the original snp file
        bed <- tempfile()
        snp %>%
          dplyr::mutate(start = pos - 1, end = pos) %>%
          dplyr::select(chrom, start, end) %>%
          dplyr::sample_n(nrow(snp) - n) %>%
          dplyr::arrange() %>%
          readr::write_tsv(bed, col_names = FALSE)
        # generate a "subset" based on thad BED file
        new_data <- filter_bed(data, bed)

        nrow(read_snp_file(new_data$exclude)) == n
    })
    expect_true(all(successes))
})

test_that("Resetting returns an EIGENSTRAT object to an original state", {
  new_data <- data
  new_data$exclude <- "exclude.snp"
  new_data$group <- "new_labels.ind"
  new_data <- reset(new_data)

  expect_true(all(sapply(c("ind", "snp", "geno"),
                         function(i) data[[i]] == new_data[[i]])))
  expect_true(is.null(new_data$groups))
  expect_true(is.null(new_data$exclude))
})
