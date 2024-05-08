#' Read an EIGENSTRAT ind/snp/geno file.
#'
#' These functions each read one part of the EIGENSTRAT dataset trio.
#'
#' Note that \code{read_geno()} will only read plain-text geno files, not compressed ones.
#'
#' @param data EIGENSTRAT data object.
#'
#' @return A data.frame object.
#'
#' @export
read_ind <- function(data) {
  path <- ifelse(is.null(data$group), data$ind, data$group)
  utils::read.table(path, col.names = c("id", "sex", "label"), stringsAsFactors = FALSE) %>%
    tibble::as_tibble(.name_repair = "minimal")
}


#' @rdname read_ind
#' @param exclude Read the list of excluded SNPs?
#'
#' @export
read_snp <- function(data, exclude = FALSE) {
  readr::read_table(
    ifelse(exclude, data$exclude, data$snp),
    col_names = c("id", "chrom", "gen", "pos", "ref", "alt"),
    col_types = "ccdicc",
    progress = FALSE
  )
}


#' @rdname read_ind
#' @export
read_geno <- function(data) {
  ind <- read_ind(data)$id

  # get the number of samples in the geno file
  n <- nchar(readLines(data$geno, 1))
  readr::read_fwf(
    data$geno,
    col_positions = readr::fwf_widths(rep(1, n), ind),
    col_types = readr::cols(.default = "i"),
    progress = FALSE
  ) %>%
    replace(., . == 9, NA)
}



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
  df <- replace(df, is.na(df), 9)
  writeLines(apply(df, 1, paste, collapse = ""), con = file)
}

