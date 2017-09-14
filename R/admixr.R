# High level wrapper functions --------------------------------------------------


#' Run an f4-ratio analysis and return the results as a data.frame.
#'
#' @param X, A, B, C, O Population names, using the terminology of
#'     Patterson et al., 2012
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix' argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @export
qpF4ratio <- function(X, A, B, C, O,
                     prefix=NULL, geno=NULL, snp=NULL, ind=NULL, badsnp=NULL,
                     dir_name=NULL) {
    check_presence(c(X, A, B, C, O), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qpF4ratio__", A, "_", B, "_", C, "_", O)
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qpF4ratio_pop_file(X=X, A=A, B=B, C=C, O=O, file=files[["pop_file"]])
    create_par_file(files[["par_file"]], files[["pop_file"]],
                    prefix, geno, snp, ind, badsnp)

    run_cmd("qpF4ratio", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qpF4ratio(files[["log_file"]]) %>% mutate(setup=setup)
}


#' Calculate D statistics or F4 statistics (which is just the
#' numerator of a D statistic) and return the results as a data.frame.
#'
#' @param W, X, Y, Z Population names, using the terminology of
#'     Patterson et al., 2012
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix' argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @export
qpDstat <- function(W, X, Y, Z,
                    prefix=NULL, geno=NULL, snp=NULL, ind=NULL, badsnp=NULL,
                    dir_name=NULL, f4mode=FALSE) {
    check_presence(c(W, X, Y, Z), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qpDstat__", W, "_", X, "_", Y, "_", Z)
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qpDstat_pop_file(W, X, Y, Z, file=files[["pop_file"]])
    create_par_file(files[["par_file"]], files[["pop_file"]],
                    prefix, geno, snp, ind, badsnp, f4mode)

    run_cmd("qpDstat", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qpF4ratio(files[["log_file"]])
}


# Reading output log files --------------------------------------------------


#' Read output log file from a qpF4ratio run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with the parsed results.
#' @export
read_qpF4ratio <- function(file) {
    log_lines <- readLines(file) %>% .[!str_detect(., "warning")]

    # extract the number of analyzed test populations/individuals
    # (corresponding to the number of rows of the results table)
    n_pops <- log_lines[which(str_detect(log_lines, "^nplist:"))] %>%
        str_extract("[0-9]+$") %>%
        as.integer

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[(length(log_lines) - n_pops) : (length(log_lines) - 1)] %>%
        str_replace("result: ", "") %>%
        str_replace_all(":", "") %>%
        str_replace_all(" +", " ") %>%
        str_replace("^ ", "")

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        read_delim(delim=" ", col_names=FALSE) %>%
        setNames(c("A", "O", "X", "C", "A", "O", "B", "C", "alpha", "stderr", "z")) %>%
        .[c("A", "B", "X", "C", "O", "alpha", "stderr", "z")]

    res_df
}


#' Read output log file from a qpDstat run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with the parsed results.
#' @export
read_qpDstat <- function(file) {
    log_lines <- readLines(file) %>% .[!str_detect(., "warning")]

    # extract the number of analyzed population quadruples
    n_quads <- length(log_lines) - (which(str_detect(log_lines, "^nrows, ncols:"))) - 1

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[(length(log_lines) - n_quads) : (length(log_lines) - 1)] %>%
        str_replace("result: ", "") %>%
        str_replace_all(" +", " ") %>%
        str_replace_all("^ | $", "")

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        read_delim(delim=" ", col_names=FALSE) %>%
        .[c(1:6, ncol(.) - 2, ncol(.) - 1, ncol(.))] %>% # remove column with "best" if present
        setNames(c("W", "X", "Y", "Z", "Dstat", "Zscore", "BABA", "ABBA", "n_snps"))

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
        lines <- str_replace(lines, regex, merge_into)
    }

    writeLines(lines, modified_file)
}


#' Read an EIGENSTRAT 'ind' file.
#'
#' @param file Path to the file.
#'
#' @return Data frame with the sample identifier, sex and label
#'     columns.
#' @export
read_ind <- function(file) {
    read_table2(file, col_names=c("id", "sex", "label"))
}


#' Read an EIGENSTRAT 'geno' file.
#'
#' @param file Path to the file.
#'
#' @return Data frame with columns containing "genotypes" of each
#'     sample (0/1/9 as defined by the EIGENSTRAT format).
#' @export
read_geno <- function(file, inds=NULL) {
    # get the number of samples in the geno file
    n <- nchar(readLines(file, 1))
    read_fwf(file, col_positions=fwf_widths(rep(1, n), inds))
}


#' Read an EIGENSTRAT 'snp' file.
#'
#' @param file Path to the file.
#'
#' @return Data frame with information about each SNP (columns defined by the EIGENSTRAT format).
#' @export
read_snp <- function(snp_file) {
    read_table(snp_file, col_names=c("id", "chrom", "gen", "pos", "alt", "ref"))

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
    summarise_all(geno, funs(fn(. != 9)))
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
    summarise_all(geno, funs(fn(. == 9)))
}

