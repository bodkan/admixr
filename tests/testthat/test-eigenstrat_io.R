context("EIGENSTRAT file I/O functionality")

prefix <- file.path(admixtools_path(), "convertf", "example")
file.copy(from = paste0(prefix, ".eigenstratgeno"), to = paste0(prefix, ".geno"), overwrite = TRUE)

# read_ind / write_ind ----------------------------------------------------

test_that("read_ind and write_ind are inverse functions", {
  true <- read_ind(paste0(prefix, ".ind"))

  new_file <- tempfile()
  write_ind(true, new_file)

  new <- read_ind(new_file)

  expect_equal(true, new)
})

# read_geno / write_geno --------------------------------------------------

test_that("read_geno and write_geno are inverse functions", {
  true <- read_geno(paste0(prefix, ".geno"))
  
  tmp <- tempfile()
  new_file <- paste0(tmp, ".geno")
  file.copy(paste0(prefix, ".ind"), paste0(tmp, ".ind"))
  write_geno(true, new_file)
  
  new <- read_geno(new_file)
  
  expect_equal(true, new)
})

# read_snp / write_snp ----------------------------------------------------

test_that("read_snp and write_snp are inverse functions", {
  true <- read_snp(paste0(prefix, ".snp"))
  
  new_file <- tempfile()
  write_snp(true, new_file)
  
  new <- read_snp(new_file)
  
  expect_equal(true, new)
})

# reading / writing EIGENSTRAT data at once -------------------------------

test_that("read_eigenstrat and write_eigenstrat are inverse functions", {
  true <- read_eigenstrat(prefix)
  
  prefix_dir <- tempdir()
  new_prefix <- file.path(prefix_dir, "test")
  write_eigenstrat(prefix = new_prefix, ind = true$ind, snp = true$snp, geno = true$geno)

  new <- read_eigenstrat(new_prefix)
  
  expect_equal(true, new)
})
