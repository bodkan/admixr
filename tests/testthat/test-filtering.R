test_that("Potential aDNA SNPs are correctly removed", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(admixtools_path(), "data/allmap"))
  output <- tempfile()
  system(
    sprintf("awk '($5 == \"C\" && $6 == \"T\") ||
                  ($5 == \"T\" && $6 == \"C\") ||
                  ($5 == \"G\" && $6 == \"A\") ||
                   ($5 == \"A\" && $6 == \"G\")' %s > %s", data$snp, output)
  )
  n_shell <- nrow(read_snp_file(output))
  n_admixr <- nrow(read_snp(transversions_only(data), exclude = TRUE))
  expect_equal(n_shell, n_admixr)
})

if (admixtools_present()) {
  orig_data <- eigenstrat(
      prefix = file.path(admixtools_path(), "convertf", "example"),
      geno = file.path(admixtools_path(), "convertf", "example.eigenstratgeno")
  )

  # ADMIXTOOLS example data is broken and it's first SNP has a position 0,
  # although snp files have to be 0-based - let's remove the first SNP entirely
  data <- orig_data
  data$snp <- tempfile()
  data$geno <- tempfile()
  system(sprintf("tail -n+2 %s > %s", orig_data$snp, data$snp))
  system(sprintf("tail -n+2 %s > %s", orig_data$geno, data$geno))
}

test_that("filter_sites correctly fails at no overlap", {
  skip_on_cran()
  skip_on_os("windows")

  # create a BED file that has the same positions as the original EIGENSTRAT
  bed <- tempfile()
  snp_to_bed(data$snp, bed)

  # verify that no overlaps leads to error
  expect_error(filter_bed(data, bed, remove = TRUE))
})

test_that("filter_sites correctly handles complete overlap", {
  skip_on_cran()
  skip_on_os("windows")

  # create a BED file that has the same positions as the original EIGENSTRAT
  bed <- tempfile()
  snp_to_bed(data$snp, bed)

  # generate a "subset" based on that BED file
  output <- tempfile()
  # suppress weird tibble deprecation warnings
  # readr still calls data_frame() function although it is deprecated now
  suppressWarnings(new_data <- filter_bed(data = data, bed = bed))

  # verify that both EIGENSTRAT datasets are the same
  expect_true(nrow(read_snp_file(new_data$exclude)) == 0)
})

test_that("Overlap returns a correct number of sites", {
  skip_on_cran()
  skip_on_os("windows")

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
    # generate a "subset" based on that BED file
    # FIX WARNINGS due to deprecation (see also above)
    suppressWarnings(new_data <- filter_bed(data = data, bed = bed))

    nrow(read_snp_file(new_data$exclude)) == n
  })
  expect_true(all(successes))
})

test_that("Resetting returns an EIGENSTRAT object to an original state", {
  skip_on_cran()
  skip_on_os("windows")

  new_data <- data
  new_data$exclude <- "exclude.snp"
  new_data$group <- "new_labels.ind"
  new_data <- reset(new_data)

  expect_true(all(sapply(c("ind", "snp", "geno"),
                         function(i) data[[i]] == new_data[[i]])))
  expect_true(is.null(new_data$groups))
  expect_true(is.null(new_data$exclude))
})

test_that("Filtering works when data is piped into a calculation", {
  skip_on_cran()
  skip_on_os("windows")

  data_dir <- file.path(admixtools_path(), "data")
  examples_dir <- file.path(admixtools_path(), "examples")

  data <- eigenstrat(file.path(data_dir, "allmap"))
  pops <- read_pops(file.path(examples_dir, "list_qpDstat1"),
                    c("W", "X", "Y", "Z"))
  expect_silent(
    data %>%
      transversions_only %>%
      d(W = pops$W, X = pops$X, Y = pops$Y, Z = pops$Z)
  )
})

test_that("filter_bed correctly handles strange contig names", {
  skip_on_cran()
  skip_on_os("windows")

  data <- eigenstrat(file.path(admixtools_path(), "data/allmap"))

  set.seed(123)

  # these are the chromosomes observed in the snp file
  # unique(read_snp(data)$chrom)

  # create a BED file that has a mixture of positions in the original snp file,
  # but also has positions with weird contigs
  # originally pointed out here as an issue which doesn't seem to be happening:
  # https://github.com/bodkan/admixr/issues/72
  bed_file <- tempfile()
  modified_snp <- tempfile()
  read_snp(data) %>%
    dplyr::sample_n(10) %>%
    dplyr::arrange(chrom, pos) %>%
    dplyr::mutate(chrom = c(chrom[1:5], paste0("asdfchr", chrom[6:10]))) %>%
    write_snp(modified_snp)
  snp_to_bed(modified_snp, bed_file)

  # check the generated BED file
  # readr::read_tsv(bed_file, col_names = c("chrom", "start", "end"), show_col_types = FALSE)

  filtered_data1 <- filter_bed(data, bed_file, remove = TRUE)
  expect_true(nrow(read_snp_file(filtered_data1$exclude)) == 5)

  # an optional argument for bedtools has been implemented, allowing `-nonamecheck` to
  # silence warnings about inconsistent naming conventions -- again, motivated by the
  # issue linked above but it seems this has no practical benefit (I will keep the
  # `bedtools_args = ""` option to allow `-sorted` to be used with large BED files)
  bed_file <- tempfile()
  modified_snp <- tempfile()
  read_snp(data) %>%
    dplyr::sample_n(10) %>%
    dplyr::mutate(chrom = c(chrom[1:5], paste0("chr", chrom[6:10]))) %>%
    dplyr::arrange(chrom, pos) %>%
    write_snp(modified_snp)
  snp_to_bed(modified_snp, bed_file)

  # again, check the produced BED file
  # readr::read_tsv(bed_file, col_names = c("chrom", "start", "end"), show_col_types = FALSE)

  # filtered_data2 <- filter_bed(data, bed_file, remove = TRUE)
  filtered_data3 <- filter_bed(data, bed_file, remove = TRUE, bedtools_args = "-nonamecheck")

  # expect_true(nrow(read_snp_file(filtered_data2$exclude)) == 5)
  expect_true(nrow(read_snp_file(filtered_data3$exclude)) == 5)
})
