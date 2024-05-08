if (admixtools_present()) {
  prefix <- file.path(admixtools_path(), "convertf", "example")
  data <- eigenstrat(prefix, geno = paste0(prefix, ".eigenstratgeno"))
}

# read_ind / write_ind ----------------------------------------------------

test_that("read_ind and write_ind are inverse functions", {
  skip_on_cran()
  skip_on_os("windows")

  true <- read_ind(data)

  new_file <- tempfile()
  write_ind(true, new_file)

  new <- readr::read_tsv(new_file, col_names = c("id", "sex", "label"))

  expect_equivalent(true, new)
})

# read_geno / write_geno --------------------------------------------------

test_that("read_geno and write_geno are inverse functions", {
  skip_on_cran()
  skip_on_os("windows")

  true <- read_geno(data)
  
  tmp <- tempfile()
  new_file <- paste0(tmp, ".geno")
  file.copy(paste0(prefix, ".ind"), paste0(tmp, ".ind"))
  write_geno(true, new_file)
  
  expect_equal(readLines(data$geno), readLines(new_file))
})

# read_snp / write_snp ----------------------------------------------------

test_that("read_snp and write_snp are inverse functions", {
  skip_on_cran()
  skip_on_os("windows")

  true <- read_snp(data)
  
  new_file <- tempfile()
  write_snp(true, new_file)
  
  new <- readr::read_table(
    new_file,
    col_types = "ccdicc",
    col_names = c("id", "chrom", "gen", "pos", "ref", "alt")
  )
  
  expect_equal(true, new)
})

