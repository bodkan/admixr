admixtools_present <- function() {
  Sys.which("qpDstat") != ""
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

admixtools_path <- function() {
    ## ugly hack to enable testing on my macOS where I have ADMIXTOOLS
    ## binaries symlinked to ~/local/bin
    if (system("uname", intern = TRUE) == "Darwin") {
        return("~/local/AdmixTools-6.0/")
    }

    stringr::str_replace(Sys.which("qpDstat"), "/bin.*", "")
}
