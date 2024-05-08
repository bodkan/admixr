# These tests run the same calculations as those performed by ADMIXTOOLS'
# built in tests (see test-01_admixtools_cmds.R), but they do so using
# admixr's wrapper functions. The results of both tests (i.e. contents of
# output log files) are then compared to each other.

if (admixtools_present()) {
    snps <- eigenstrat(download_data())
}

test_that("qpAdm_rotation() `fulloutput = TRUE` object is of the right format", {
    skip_on_cran()
    skip_on_os("windows")

    models1 <- qpAdm_rotation(
        data = snps,
        target = "French",
        candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai"),
        minimize = FALSE,
        nsources = 2,
        ncores = 1,
        fulloutput = TRUE
    )
    expect_equal(sort(names(models1)), sort(c("proportions", "ranks", "subsets")))
})


test_that("qpAdm_rotation() `fulloutput = FALSE` object is of the right format", {
    skip_on_cran()
    skip_on_os("windows")

    models2 <- qpAdm_rotation(
        data = snps,
        target = "French",
        candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai"),
        minimize = FALSE,
        nsources = 2,
        ncores = 1,
        fulloutput = FALSE
    )
    expect_true(is.data.frame(models2))
})


test_that("qpAdm_rotation() output object carries all required annotation", {
    skip_on_cran()
    skip_on_os("windows")

    models1 <- qpAdm_rotation(
        data = snps,
        target = "French",
        candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai"),
        minimize = FALSE,
        nsources = 2,
        ncores = 1,
        fulloutput = TRUE
    )
    expect_true(!is.null(attr(models1, "log_output")))
})


test_that("qpAdm_filter() filters p-value and proportions", {
    skip_on_cran()
    skip_on_os("windows")
    models1 <- qpAdm_rotation(
        data = snps,
        target = "French",
        candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai"),
        minimize = FALSE,
        nsources = 2,
        ncores = 1,
        fulloutput = TRUE
    )
    filtered <- qpAdm_filter(models1, p = 0.01)
    expect_true(all(c(filtered$proportions$pvalue > 0.01,
                      filtered$proportions$prop1 < 1,
                      filtered$proportions$prop1 > 0,
                      filtered$proportions$prop2 < 1,
                      filtered$proportions$prop2 > 0)))
})

test_that("qpAdm_fiter works only on objects of qpAdm_rotation type", {
    x <- data.frame()
    class(x) <- "random"
    expect_error(qpAdm_filter(x))
})

test_that("qpAdm gives the correct number of source/stderr columns (n = 2)", {
  skip_on_cran()
  skip_on_os("windows")

  # two sources
  nsources <- 2
  sources2 <- qpAdm_rotation(
    data = snps,
    target = "French",
    candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai"),
    minimize = FALSE,
    nsources = nsources,
    ncores = 1,
    fulloutput = TRUE
  )
  expect_true(all(unique(sources2$ranks$rank) == (nsources - 1) : nsources))
  expect_true(all(grep("source", colnames(sources2$proportions), value = TRUE) == paste0("source", 1:nsources)))

  # three sources
  nsources <- 3
  sources3 <- qpAdm_rotation(
    data = snps,
    target = "French",
    candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai", "Denisova", "Sardinian"),
    minimize = FALSE,
    nsources = nsources,
    ncores = 1,
    fulloutput = TRUE
  )
  expect_true(all(unique(sources3$ranks$rank) == (nsources - 1) : nsources))
  expect_true(all(grep("source", colnames(sources3$proportions), value = TRUE) == paste0("source", 1:nsources)))

  # four sources
  nsources <- 4
  sources4 <- qpAdm_rotation(
    data = snps,
    target = "French",
    candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija", "Altai", "Sardinian", "Khomani_San", "Chimp", "Han"),
    minimize = FALSE,
    nsources = nsources,
    ncores = 1,
    fulloutput = TRUE
  )
  expect_true(all(unique(sources4$ranks$rank) == (nsources - 1) : nsources))
  expect_true(all(grep("source", colnames(sources4$proportions), value = TRUE) == paste0("source", 1:nsources)))
})