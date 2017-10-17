# Reading EIGENSTRAT files --------------------------------------------------


#' Read an EIGENSTRAT 'ind' file.
#'
#' @param file Path to the file.
#'
#' @return Data frame with the sample identifier, sex and label
#'     columns (columns defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
read_ind <- function(file) {
    read_table2(file, col_names=c("id", "sex", "label"), col_types="ccc")
}


#' Read an EIGENSTRAT 'snp' file.
#'
#' @param file Path to the file.
#'
#' @return Data frame with information about each SNP (columns defined
#'     by the EIGENSTRAT format).
#'
#' @export
#' @import readr
read_snp <- function(snp_file) {
    read_table2(snp_file, col_names=c("id", "chrom", "gen", "pos", "ref", "alt"),
                col_types="ccdicc", progress=FALSE)
}


#' Read an EIGENSTRAT 'geno' file.
#'
#' @param file Path to the geno file.
#' @param ind_file Path to the ind file to read sample names from.
#'
#' @return Data frame with columns containing "genotypes" of each
#'     sample (0/1/9 as defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
read_geno <- function(file, ind_file=NULL) {
    if (!is.null(ind_file)) {
        inds <- read_ind(ind_file)$id
    } else {
        inds <- NULL
    }

    # get the number of samples in the geno file
    n <- nchar(readLines(file, 1))
    read_fwf(file, col_positions=fwf_widths(rep(1, n), inds),
             col_types=cols(.default="i"), progress=FALSE)
}


#' Read a tripplet of EIGENSTRAT (geno/snp/ind files) files.
#'
#' @param file Path to the file.
#'
#' @return List of three data frames (one element for geno/snp/ind).
#'
#' @export
read_eigenstrat <- function(prefix=NULL) {
    list(
        geno=read_geno(paste0(prefix, ".geno"), paste0(prefix, ".ind")),
        snp=read_snp(paste0(prefix, ".snp")),
        ind=read_ind(paste0(prefix, ".ind"))
    )
}



# Writing EIGENSTRAT files --------------------------------------------------


#' Write an EIGENSTRAT 'ind' file.
#'
#' @param file Path to the file.
#' @param ind Data frame with the sample identifier, sex and label
#'     columns (columns defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
write_ind <- function(ind_file, df) {
    write_tsv(df, ind_file, col_names=FALSE)
}


#' Write an EIGENSTRAT 'snp' file.
#'
#' @param file Path to the file.
#' @param snp Data frame with information about each SNP (columns
#'     defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
write_snp <- function(snp_file, df) {
    write_tsv(df, snp_file, col_names=FALSE)
}


#' Write an EIGENSTRAT 'geno' file.
#'
#' @param file Path to the file.
#' @param geno Data frame with columns containing "genotypes" of each
#'     sample (0/1/9 as defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
write_geno <- function(geno_file, df) {
    writeLines(apply(df, 1, paste, collapse=""), con=geno_file)
}


#' Write a tripplet of EIGENSTRAT (geno/snp/ind files) files.
#'
#' @param prefix Prefix of the geno/snp/ind files (can
#'     include the path).
#' @param ind data.frame with data in a 'ind' format
#' @param snp data.frame with data in a 'snp' format
#' @param geno data.frame with data in a 'geno' format
#'
#' @return List of three data frames (one element for geno/snp/ind).
#'
#' @export
write_eigenstrat <- function(prefix, ind, snp, geno) {
    write_ind(paste0(prefix, ".ind"), ind)
    write_snp(paste0(prefix, ".snp"), snp)
    write_geno(paste0(prefix, ".geno"), geno)
}