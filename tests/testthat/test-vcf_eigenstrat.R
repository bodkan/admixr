context("EIGENSTRAT-VCF conversion")

prefix <- file.path(admixtools_path(), "convertf", "example")
file.copy(from = paste0(prefix, ".eigenstratgeno"), to = paste0(prefix, ".geno"), overwrite = TRUE)

# read_ind / write_ind ----------------------------------------------------

test_that("vcf_to_eigenstrat and eigenstrat_to_vcf are inverse functions", {
  # create VCF from EIGENSTRAT
  new_vcf <- tempfile()
  eigenstrat_to_vcf(prefix, new_vcf, compress = FALSE, index = FALSE)

  # convert that VCF back to EIGENSTRAT
  new_prefix <- tempfile() 
  vcf_to_eigenstrat(new_vcf, new_prefix)

  # load both old and new EIGENSTRAT datasets and compare them
  # (because VCF cannot store sex, pop. labels, genetic distance, and other
  # annotation data, we will compare only the GT 0/1/2/9 tables)
  true <- read_eigenstrat(prefix)
  new <- read_eigenstrat(new_prefix)
  expect_equal(true$geno, new$geno)
})
