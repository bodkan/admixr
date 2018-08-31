# Reading EIGENSTRAT files --------------------------------------------------


#' Read an EIGENSTRAT ind/snp/geno file.
#'
#' @param file Path to a file.
#'
#' @return A data.frame object.
#'
#' @export
read_ind <- function(file) {
    readr::read_table2(file, col_names = c("id", "sex", "label"), col_types = "ccc")
}


#' @rdname read_ind
#' @export
read_snp <- function(file) {
    readr::read_table2(
      file,
      col_names = c("id", "chrom", "gen", "pos", "ref", "alt"),
      col_types = "ccdicc",
      progress = FALSE
    )
}


#' @rdname read_ind
#' @export
read_geno <- function(file) {
    ind_file = stringr::str_replace(file, "geno$", "ind")
    inds <- read_ind(ind_file)$id

    # get the number of samples in the geno file
    n <- nchar(readLines(file, 1))
    readr::read_fwf(
      file,
      col_positions = readr::fwf_widths(rep(1, n), inds),
      col_types = readr::cols(.default = "i"),
      progress = FALSE
    )
}


#' Read a trio of EIGENSTRAT geno/snp/ind files.
#'
#' @param prefix EIGENSTRAT geno/snp/ind prefix.
#'
#' @return List of three data.frame objects.
#'
#' @export
read_eigenstrat <- function(prefix = NULL) {
    list(
        geno = read_geno(paste0(prefix, ".geno")),
        snp = read_snp(paste0(prefix, ".snp")),
        ind = read_ind(paste0(prefix, ".ind"))
    )
}



# Writing EIGENSTRAT files --------------------------------------------------


#' Write an EIGENSTRAT ind/snp/geno file.
#'
#' @param df A data.frame object.
#' @param file Path to an output file.
#'
#' @export
write_ind <- function(df, file) {
    readr::write_tsv(df, file, col_names = FALSE)
}


#' @rdname write_ind
#' @export
write_snp <- function(df, file) {
    readr::write_tsv(df, file, col_names = FALSE)
}


#' @rdname write_ind
#' @export
write_geno <- function(df, file) {
    writeLines(apply(df, 1, paste, collapse = ""), con = file)
}


#' Write a trio of EIGENSTRAT geno/snp/ind files.
#'
#' @param prefix EIGENSTRAT geno/snp/ind prefix.
#' @param ind,snp,geno EIGENSTRAT data as data.frames.
#'
#' @export
write_eigenstrat <- function(prefix, ind, snp, geno) {
    write_ind(ind, paste0(prefix, ".ind"))
    write_snp(snp, paste0(prefix, ".snp"))
    write_geno(geno, paste0(prefix, ".geno"))
}
