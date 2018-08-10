# Reading EIGENSTRAT files --------------------------------------------------


#' Read an EIGENSTRAT 'ind' file.
#'
#' @param file Path to the file.
#'
#' @return A data.frame object with sample identifier, sex and label columns
#'     (as defined by the EIGENSTRAT format).
#'
#' @export
read_ind <- function(file) {
    readr::read_table2(file, col_names = c("id", "sex", "label"), col_types = "ccc")
}


#' Read an EIGENSTRAT 'snp' file.
#'
#' @param file Path to the file.
#'
#' @return A data.frame object with SNP data (columns defined by the
#'     EIGENSTRAT format).
#'
#' @export
read_snp <- function(file) {
    readr::read_table2(
      file,
      col_names = c("id", "chrom", "gen", "pos", "ref", "alt"),
      col_types = "ccdicc",
      progress = FALSE
    )
}


#' Read an EIGENSTRAT 'geno' file.
#'
#' @param file Path to the geno file.
#' @param ind_file Path to the ind file to read sample names from.
#'
#' @return A data.frame object with genotypes of each sample (0/1/9 as defined
#'     by the EIGENSTRAT format).
#'
#' @export
read_geno <- function(file, ind_file = NULL) {
    if (!is.null(ind_file)) {
        inds <- read_ind(ind_file)$id
    } else {
        inds <- NULL
    }

    # get the number of samples in the geno file
    n <- nchar(readLines(file, 1))
    readr::read_fwf(
      file,
      col_positions = readr::fwf_widths(rep(1, n), inds),
      col_types = readr::cols(.default = "i"),
      progress = FALSE
    )
}


#' Read a tripplet of EIGENSTRAT (geno/snp/ind files) files.
#'
#' @param prefix EIGENSTRAT prefix of geno/snp/ind files.
#' @param ind_suffix,snp_suffix,geno_suffix Alternative EIGENSTRAT suffixes.
#'
#' @return List of three data.frame objects (one for geno/snp/ind data).
#'
#' @export
read_eigenstrat <- function(prefix = NULL,
                            geno_suffix = ".geno",
                            ind_suffix = ".ind",
                            snp_suffix = ".snp") {
    list(
        geno = read_geno(paste0(prefix, geno_suffix), paste0(prefix, ind_suffix)),
        snp = read_snp(paste0(prefix, snp_suffix)),
        ind = read_ind(paste0(prefix, ind_suffix))
    )
}



# Writing EIGENSTRAT files --------------------------------------------------


#' Write an EIGENSTRAT 'ind' file.
#'
#' @param df A data.frame object with genotypes of each sample (0/1/9 as defined
#'     by the EIGENSTRAT format).
#' @param file Path to the file.
#'
#' @export
write_ind <- function(df, file) {
    readr::write_tsv(df, file, col_names = FALSE)
}


#' Write an EIGENSTRAT 'snp' file.
#'
#' @param df A data.frame object with SNP data (columns defined by the
#'     EIGENSTRAT format).
#' @param file Path to the file.
#'
#' @export
write_snp <- function(df, file) {
    readr::write_tsv(df, file, col_names = FALSE)
}


#' Write an EIGENSTRAT 'geno' file.
#'
#' @param df Data frame with columns containing "genotypes" of each
#'     sample (0/1/9 as defined by the EIGENSTRAT format).
#' @param file Path to the file.
#'
#' @export
write_geno <- function(df, file) {
    writeLines(apply(df, 1, paste, collapse = ""), con = file)
}


#' Write a triplet of EIGENSTRAT (geno/snp/ind files) files.
#'
#' @param prefix Prefix of the geno/snp/ind files (can
#'     include the path).
#' @param ind data.frame object with data in a 'ind' format
#' @param snp data.frame object with data in a 'snp' format
#' @param geno data.frame object with data in a 'geno' format
#'
#' @export
write_eigenstrat <- function(prefix, ind, snp, geno) {
    write_ind(ind, paste0(prefix, ".ind"))
    write_snp(snp, paste0(prefix, ".snp"))
    write_geno(geno, paste0(prefix, ".geno"))
}
