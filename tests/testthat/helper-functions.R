admixtools_present <- function() {
  as.vector(Sys.which("qpDstat") != "")
}

read_pops <- function(filename, columns) {
  df <- setNames(read.table(filename, stringsAsFactors = FALSE), columns)
  lapply(df, function(col) unique(col))
}

read_snp_file <- function(path) {
  readr::read_table(
    path,
    col_types = "ccdicc",
    col_names = c("id", "chrom", "gen", "pos", "ref", "alt")
  )
}

snp_to_bed <- function(snp, bed) {
  read_snp_file(snp) %>%
        dplyr::mutate(start = pos - 1, end = pos) %>%
        dplyr::select(chrom, start, end) %>%
        readr::write_tsv(bed, col_names = FALSE)
}

# Find the path to the root of the ADMIXTOOLS directory
admixtools_path <- function() {
  path <- Sys.which("qpDstat")
  is_symlink <- base::isTRUE(base::nzchar(base::Sys.readlink(path), keepNA = TRUE))
  if (is_symlink) path <- Sys.readlink(path)
  # if ADMIXTOOLS was set up with `make install`, the binary will be under
  # ./bin, otherwise it will be under ./src -- either way, remove those prefixes
  # to get the root of the ADMIXTOOLS source directory
  stringr::str_replace(path, "(/src\\/qpDstat|/bin\\/qpDstat)", "")
}
