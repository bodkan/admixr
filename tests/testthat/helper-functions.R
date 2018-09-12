admixtools_present <- function() {
  system("which qpDstat", ignore.stdout = TRUE) == 0
}

read_pops <- function(filename, columns) {
  df <- setNames(read.table(filename, stringsAsFactors = FALSE), columns)
  lapply(df, function(col) unique(col))
}

read_snp_file <- function(path) {
  readr::read_table2(
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

