#' Convert EIGENSTRAT files into a VCF file.
#'
#' This function reads the genotypes from a 'geno' file, their
#' coordinates and reference/alternative alleles and the individuals'
#' names from an 'ind' file and generates a minimalistic VCF file.
#'
#' Compressing and indexing requires having "bgzip" and "tabix"
#' commands in $PATH.
#'
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param vcf_file Path to the VCF file that will be generated.
#' @param compress Compress the VCF with bgzip?
#' @param index Index the VCF with tabix?
#'
#' @export
eigenstrat_to_vcf <- function(prefix, vcf_file, compress=TRUE, index=TRUE) {
    geno <- read_geno(paste0(prefix, ".geno"), paste0(prefix, ".ind"))
    snp <- read_snp(paste0(prefix, ".snp"))
    ind <- read_ind(paste0(prefix, ".ind"))

    # construct a minimal VCF header
    header <- c("##fileformat=VCFv4.1",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                sapply(unique(snp$chrom), function(chrom) {
                    paste0("##contig=<ID=", chrom, ">")
                }))

    # generate a dataframe with the "body" of the VCF file (info and GT columns)
    info_cols <-
      dplyr::mutate(snp, ID=".", QUAL="0", FILTER=".", INFO=".", FORMAT="GT") %>%
      dplyr::select(`#CHROM`=chrom, POS=pos, ID, REF=ref, ALT=alt, QUAL, FILTER, INFO, FORMAT)
    gt_cols <- mutate_all(geno, eigenstrat_to_gt)
    body_cols <- bind_cols(info_cols, gt_cols)

    writeLines(header, vcf_file)
    readr::write_tsv(body_cols, vcf_file, col_names=TRUE, append=TRUE)

    if (compress) system(paste("bgzip", vcf_file))
    if (index) system(paste("tabix", paste0(vcf_file, ".gz")))
}


#' Convert VCF file into three-file EIGENSTRAT format.
#'
#' @param vcf_file Path to the VCF file.
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path) that will be generated.
#'
#' @export
vcf_to_eigenstrat <- function(vcf_file, prefix) {
    vcf <- readr::read_tsv(vcf_file, comment="##") %>%
        dplyr::rename(chrom=`#CHROM`, pos=POS, ref=REF, alt=ALT) %>%
        dplyr::select(-c(ID, QUAL, FILTER, INFO, FORMAT))

    # generate dataframes with the 3 EIGENSTRAT info tables
    snp <- dplyr::select(vcf, chrom, pos, ref, alt) %>%
        dplyr::mutate(snp_id=paste(chrom, pos, sep="_"), gen_dist="0.0") %>%
        dplyr::select(snp_id, chrom, gen_dist, pos, ref, alt)
    ind <- tibble::tibble(
        sample_id=select(vcf, -c(chrom, pos, ref, alt)) %>% names,
        sex="U",
        label=sample_id
    )
    geno <- dplyr::select(vcf, -c(chrom, pos, ref, alt)) %>%
        dplyr::mutate_all(gt_to_eigenstrat)

    # write all three EIGENSTRAT files
    readr::write_tsv(snp, paste0(prefix, ".snp"), col_names=FALSE)
    readr::write_tsv(ind, paste0(prefix, ".ind"), col_names=FALSE)
    writeLines(apply(geno, 1, paste, collapse=""), paste0(prefix, ".geno"))
}


# Convert VCF-like GT string(s) into EIGENSTRAT genotypes.
# gt_to_eigenstrat(c(".|.", "./.", ".", "0|0", "0/0", "0", "0|1", "1|0", "0/1", "1|1", "1/1", "1"))
gt_to_eigenstrat <- function(gts) {
    eigen_gts <- gts %>%
        stringr::str_replace("^0$",   "2") %>%
        stringr::str_replace("^1$",   "0") %>%
        stringr::str_replace("^\\.\\|\\.$", "9") %>%
        stringr::str_replace("^\\./\\.$", "9") %>%
        stringr::str_replace("^\\.$",   "9") %>%
        stringr::str_replace("^0\\|0$", "2") %>%
        stringr::str_replace("^0/0$", "2") %>%
        stringr::str_replace("^0\\|1$", "1") %>%
        stringr::str_replace("^1\\|0$", "1") %>%
        stringr::str_replace("^0/1$", "1") %>%
        stringr::str_replace("^1\\|1$", "0") %>%
        stringr::str_replace("^1/1$", "0")

    eigen_gts
}

# Convert VCF-like GT string(s) into EIGENSTRAT genotypes.
# eigenstrat_to_gt(c(0, 1, 2, 9, 9, 2, 1, 0))
eigenstrat_to_gt <- function(eigenstrat_gts) {
    vcf_gts <- sapply(eigenstrat_gts, function(i) {
        if      (i == 0) "1/1"
        else if (i == 1) "0/1"
        else if (i == 2) "0/0"
        else             "./."
    })

    vcf_gts
}
