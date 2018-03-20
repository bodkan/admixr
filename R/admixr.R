# Reading output log files --------------------------------------------------


#' Read output log file from a qpF4ratio run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with the parsed results.
#' @export
read_qpF4ratio <- function(file) {
    log_lines <- readLines(file) %>% .[!stringr::str_detect(., "warning")]

    # extract the number of analyzed test populations/individuals
    # (corresponding to the number of rows of the results table)
    n_pops <- log_lines[which(stringr::str_detect(log_lines, "^nplist:"))] %>%
        stringr::str_extract("[0-9]+$") %>%
        as.integer

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[(length(log_lines) - n_pops) : (length(log_lines) - 1)] %>%
        stringr::str_replace("result: ", "") %>%
        stringr::str_replace_all(":", "") %>%
        stringr::str_replace_all(" +", " ") %>%
        stringr::str_replace("^ ", "")

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        readr::read_delim(delim=" ", col_names=FALSE) %>%
        setNames(c("A", "O", "X", "C", "A", "O", "B", "C", "alpha", "stderr", "Zscore")) %>%
        .[c("A", "B", "X", "C", "O", "alpha", "stderr", "Zscore")]

    res_df
}


#' Read output log file from a qpDstat run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with the parsed results.
#' @export
read_qpDstat <- function(file) {
    log_lines <- readLines(file) %>% .[!stringr::str_detect(., "warning")]

    # extract the number of analyzed population quadruples
    n_quads <- length(log_lines) - (which(stringr::str_detect(log_lines, "^nrows, ncols:"))) - 1

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[(length(log_lines) - n_quads) : (length(log_lines) - 1)] %>%
        stringr::str_replace("result: ", "") %>%
        stringr::str_replace_all(" +", " ") %>%
        stringr::str_replace_all("^ | $", "")

    result_col <- ifelse(any(stringr::str_detect(log_lines, "f4mode: YES")), "f4", "D")

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        readr::read_delim(delim=" ", col_names=FALSE) %>%
        .[c(1:6, ncol(.) - 2, ncol(.) - 1, ncol(.))] %>% # remove column with "best" if present
        setNames(c("W", "X", "Y", "Z", result_col, "Zscore", "BABA", "ABBA", "n_snps"))

    res_df
}


#' Read output log file from a qp3Pop run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with parsed results.
#' @export
read_qp3Pop <- function(file) {
    log_lines <- readLines(file)

    # extract the number of analyzed population quadruples
    n_quads <- length(log_lines) - (which(stringr::str_detect(log_lines, "^nrows, ncols:"))) - 1

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[stringr::str_detect(log_lines, "result:")] %>%
        stringr::str_replace("result: ", "") %>%
        stringr::str_replace_all(" +", " ") %>%
        stringr::str_replace_all("^ | $", "")

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        readr::read_delim(delim=" ", col_names=FALSE) %>%
        setNames(c("A", "B", "C", "f3", "stderr", "Zscore", "n_snps"))

    res_df
}


# EIGENSTRAT manipulation utilities --------------------------------------------------


#' Merge populations from an EIGENSTRAT "ind" file under a single
#' population label.
#'
#' @param file EIGENSTRAT ind file to modify.
#' @param modified_file Modified EIGENSTRAT ind filename.
#' @param merge List of labels to merge. List names specified labels
#'     to merge into.
#' @export
merge_pops <- function(file, modified_file, merge) {
    # merge=list(ancient_NearEast=merge_what, present_NearEast=c("Yemenite_Jew", "Jordan", "Samaritan", "Bedouin", "Palestinian"))
    lines <- readLines(file)

    # iterate over the lines in the "ind" file, replacing population
    # labels with their substitutes
    for (merge_into in names(merge)) {
        regex <- paste0("(", paste(merge[[merge_into]], collapse="|"), ")$")
        lines <- stringr::str_replace(lines, regex, merge_into)
    }

    writeLines(lines, modified_file)
}



# Filtering functions  --------------------------------------------------


#' Calculate the number (or proportion) of sites with an allele
#' present (i.e. not 9) for each sample.
#'
#' @param geno EIGENSTRAT geno dataframe.
#' @param prop Calculate the proportion of non-missing alleles
#'     instead.
#'
#' @return A named vector of counts or proportions.
#' @export
snps_present <- function(geno, prop=FALSE) {
    fn <- ifelse(prop, mean, sum)
    dplyr::summarise_all(geno, funs(fn(. != 9)))
}


#' Calculate the number (or proportion) of sites with an allele
#' missing for each sample.
#'
#' @param geno EIGENSTRAT geno dataframe.
#' @param prop Calculate the proportion of missing alleles
#'     instead.
#'
#' @return A named vector of counts or proportions.
#' @export
snps_missing <- function(geno, prop=FALSE) {
    fn <- ifelse(prop, mean, sum)
    dplyr::summarise_all(geno, funs(fn(. == 9)))
}



#' Create a new set of EIGENSTRAT files by intersecting the original
#' data with a given set of coordinates.
#'
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param out_prefix Prefix of the generated EIGENSTRAT files with the
#'     subset of the data.
#' @param bed_file Path to the 3 column BED file to intersect with.
#' @param complement Perform an intersect or a complement operation?
#'
#' @export
subset_sites <- function(prefix, out_prefix, bed_file, complement=FALSE) {
    coords <- readr::read_table2(
      bed_file,
      col_names=c("chrom", "start", "pos"),
      col_types="cii",
      progress=FALSE
    ) %>%
      dplyr::select(-start)

    geno <- read_geno(paste0(prefix, ".geno"))
    snp <- read_snp(paste0(prefix, ".snp"))
    combined <- dplyr::bind_cols(snp, geno)

    # determine which function to call on the coordinates
    fun <- ifelse(complement, dplyr::anti_join, dplyr::inner_join)
    combined_subset <- fun(combined, coords, by=c("chrom", "pos"))

    # write the new snp file
    dplyr::select(combined_subset, id:alt) %>%  
      readr::write_tsv(path=paste0(out_prefix, ".snp"), col_names=FALSE)
    # write the new geno file
    dplyr::select(combined_subset, -(id:alt)) %>% 
      apply(1, paste, collapse="") %>%
      writeLines(con=paste0(out_prefix, ".geno"))
    # write the new ind file
    invisible(file.copy(from=paste0(prefix, ".ind"),
                        to=paste0(out_prefix, ".ind")))
}

