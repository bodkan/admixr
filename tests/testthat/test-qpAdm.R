# These tests run the same calculations as those performed by ADMIXTOOLS'
# built in tests (see test-01_admixtools_cmds.R), but they do so using
# admixr's wrapper functions. The results of both tests (i.e. contents of
# output log files) are then compared to each other.

context("qpAdm rotation functionality")

if (admixtools_present()) {
    snps <- eigenstrat(download_data())
    ## modeling with full output
    models1 <- qpAdm_rotation(
        data = snps,
        target = "French",
        candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai"),
        minimize = FALSE,
        nsources = 2,
        ncores = 1,
        fulloutput = TRUE
    )
    ## modelling with minimal output
    models2 <- qpAdm_rotation(
        data = snps,
        target = "French",
        candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai"),
        minimize = FALSE,
        nsources = 2,
        ncores = 1,
        fulloutput = FALSE
    )
}

test_that("qpAdm_rotation() `fulloutput = TRUE` object is of the right format", {
    skip_on_cran()
    skip_on_os("windows")
    expect_equal(sort(names(models1)), sort(c("proportions", "ranks", "subsets")))
})


test_that("qpAdm_rotation() `fulloutput = FALSE` object is of the right format", {
    skip_on_cran()
    skip_on_os("windows")
    expect_true(is.data.frame(models2))
})


test_that("qpAdm_rotation() output object carries all required annotation", {
    skip_on_cran()
    skip_on_os("windows")
    expect_true(!is.null(attr(models1, "log_output")))
})


test_that("qpAdm_filter() filters p-value and proportions", {
    skip_on_cran()
    skip_on_os("windows")
    filtered <- qpAdm_filter(models1, p = 0.01)
    expect_true(all(c(filtered$proportions$pvalue > 0.01,
                      filtered$proportions$prop1 < 1,
                      filtered$proportions$prop1 > 0,
                      filtered$proportions$prop2 < 1,
                      filtered$proportions$prop2 > 0)))
})


test_that("It is possible to switch off outgroup number minimization", {
    models <- qpAdm_rotation(
        data = snps,
        target = "French",
        candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai"),
        minimize = FALSE,
        nsources = 2,
        ncores = 1,
        fulloutput = TRUE
    )
    expect_true(unique(models$proportions$noutgroups) == 3)
})

