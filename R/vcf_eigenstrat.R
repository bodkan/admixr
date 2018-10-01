#' Conversion between VCF and EIGENSTRAT file formats
#'
#' These functions parse input in either VCF or EIGENSTRAT format and convert it.
#'
#' Please note that processing of large data sets (such as whole-genome data)
#' using these functions will be very inefficient, because they first load all
#' data into memory first, perform the conversion, and only then save the
#' output file.  
#'
#' Compressing and indexing is performed by default, and requires having
#' bgzip and tabix commands in $PATH.
#'
#' @param vcf Path to a VCF file.
#' @param eigenstrat EIGENSTRAT data object.
#' @param compress Compress the VCF with bgzip?
#' @param index Index the VCF with tabix?
#'
#' @export
xvcf_to_eigenstrat <- function(vcf, eigenstrat) {
  vcf <- readr::read_tsv(vcf, comment = "##", progress = FALSE) %>%
    dplyr::rename(chrom = `#CHROM`, pos = POS, snp_id = ID, ref = REF, alt = ALT) %>%
    dplyr::select(-c(QUAL, FILTER, INFO, FORMAT))

  # generate dataframes with the 3 EIGENSTRAT info tables
  snp <- dplyr::select(vcf, chrom, pos, snp_id, ref, alt) %>%
    dplyr::mutate(gen_dist = "0.0") %>%
    dplyr::select(snp_id, chrom, gen_dist, pos, ref, alt)
  ind <- tibble::tibble(
    sample_id = dplyr::select(vcf, -c(chrom, pos, snp_id, ref, alt)) %>% names,
    sex = "U",
    label = sample_id
  )
  geno <- dplyr::select(vcf, -c(chrom, pos, snp_id, ref, alt)) %>%
    dplyr::mutate_all(gt_to_eigenstrat)

  # write all three EIGENSTRAT files
  readr::write_tsv(snp, eigenstrat$snp, col_names = FALSE)
  readr::write_tsv(ind, eigenstrat$ind, col_names = FALSE)
  writeLines(apply(geno, 1, paste, collapse = ""), eigenstrat$geno)
}



#' @rdname vcf_to_eigenstrat
#' @export
xeigenstrat_to_vcf <- function(eigenstrat, vcf, compress = TRUE, index = TRUE) {
  geno <- read_geno(eigenstrat)
  snp <- read_snp(eigenstrat)

  # construct a minimal VCF header
  header <- c(
    "##fileformat=VCFv4.1",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    sapply(unique(snp$chrom), function(chrom) { paste0("##contig=<ID=", chrom, ">") })
  )

  # generate a dataframe with the "body" of the VCF file (info and GT columns)
  info_cols <-
    dplyr::mutate(snp, ID = snp$id, QUAL = "0", FILTER = ".", INFO = ".", FORMAT = "GT") %>%
    dplyr::select(`#CHROM` = chrom, POS = pos, ID, REF = ref, ALT = alt, QUAL, FILTER, INFO, FORMAT)
  gt_cols <- dplyr::mutate_all(geno, eigenstrat_to_gt)
  body_cols <- dplyr::bind_cols(info_cols, gt_cols)

  writeLines(header, vcf)
  readr::write_tsv(body_cols, vcf, col_names = TRUE, append = TRUE)

  if (compress) system(paste("bgzip", vcf))
  if (index) system(paste("tabix", paste0(vcf, ".gz")))
}



# Convert VCF-like GT string(s) into EIGENSTRAT genotypes.
# gt_to_eigenstrat(c(".|.", "./.", ".", "0|0", "0/0", "0", "0|1", "1|0", "0/1", "1|1", "1/1", "1"))
gt_to_eigenstrat <- function(gts) {
  gts <- stringr::str_replace(gts, "^0$",         "2")
  gts <- stringr::str_replace(gts, "^1$",         "0")
  gts <- stringr::str_replace(gts, "^\\.\\|\\.$", "9")
  gts <- stringr::str_replace(gts, "^\\./\\.$",   "9")
  gts <- stringr::str_replace(gts, "^\\.$",       "9")
  gts <- stringr::str_replace(gts, "^0\\|0$",     "2")
  gts <- stringr::str_replace(gts, "^0/0$",       "2")
  gts <- stringr::str_replace(gts, "^0\\|1$",     "1")
  gts <- stringr::str_replace(gts, "^1\\|0$",     "1")
  gts <- stringr::str_replace(gts, "^0/1$",       "1")
  gts <- stringr::str_replace(gts, "^1\\|1$",     "0")
  gts <- stringr::str_replace(gts, "^1/1$",       "0")

  gts
}



# Convert VCF-like GT string(s) into EIGENSTRAT genotypes.
eigenstrat_to_gt <- function(eigenstrat_gts) {
  vcf_gts <- sapply(eigenstrat_gts, function(i) {
    if      (i == 0) "1/1"
    else if (i == 1) "0/1"
    else if (i == 2) "0/0"
    else             "./."
  })

  vcf_gts
}
