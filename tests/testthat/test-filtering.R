context("Filtering of sites")

prefix <- file.path(admixtools_path(), "convertf", "example")
snp <- paste0(prefix, ".snp")

write_bed <- function(snp, bed) {
    read_snp(snp) %>%
        dplyr::mutate(start = pos - 1, end = pos) %>%
        dplyr::select(chrom, start, end) %>%
        readr::write_tsv(bed, col_names = FALSE)
}

test_that("filter_sites correctly handles complete overlap", {
    # create a BED file that has the same positions as the original EIGENSTRAT
    bed <- tempfile()
    write_bed(snp, bed)

    # generate a "subset" based on that BED file
    output <- tempfile()
    filter_sites(prefix, bed, output)

    # verify that both EIGENSTRAT datasets are the same
    orig_snp <- read_snp(snp)
    output_snp <- read_snp(output)
    expect_equal(orig_snp, output_snp)
})

test_that("filter_sites correctly fails at no overlap", {
    # create a BED file that has the same positions as the original EIGENSTRAT
    bed <- tempfile()
    write_bed(snp, bed)

    # verify that no overlaps leads to error
    expect_error(filter_sites(prefix, bed, "blah", remove = TRUE))
})

test_that("Overlap returns a correct number of sites", {
    snp <- read_snp(snp)

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
        filter_sites(prefix, bed, output)
        output_snp <- read_snp(output)

        nrow(output_snp) == n
    })
    expect_true(all(successes))
})

