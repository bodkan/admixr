context("EIGENSTRAT-VCF conversion")

test_that("vcf_to_eigenstrat and eigenstrat_to_vcf are inverse functions", {
  data <- eigenstrat(prefix = file.path(admixtools_path(), "convertf", "example"),
                     geno = file.path(admixtools_path(), "convertf", "example.eigenstratgeno"))
  
  # create VCF from EIGENSTRAT
  new_vcf <- tempfile()
  eigenstrat_to_vcf(eigenstrat = data, vcf = new_vcf, compress = FALSE, index = FALSE)

  # convert that VCF back to EIGENSTRAT
  new_data <- data
  new_data$ind <- tempfile()
  new_data$snp <- tempfile()
  new_data$geno <- tempfile();

  vcf_to_eigenstrat(new_vcf, new_data)

  # load both old and new EIGENSTRAT datasets and compare them
  # (because VCF cannot store sex, pop. labels, genetic distance, and other
  # annotation data, we will compare only the GT 0/1/2/9 tables)

  true <- read_geno(data)
  new <- read_geno(new_data)
  expect_equal(true, new)
})
