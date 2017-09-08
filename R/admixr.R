# High level wrapper functions ====================================================


#' Run an f4-ratio analysis and return the results as a data.frame.
#'
#' @param X, A, B, C, O Population names, using the terminology of
#'     Patterson et al., 2012
#' @param eigenstrat_prefix Prefix of the geno/snp/ind files (can
#'     include the path). If specified, geno_file/snp_file/ind_file
#'     have to be NULL and vice versa.
#' @param geno_file Path to the genotype file.
#' @param snp_file Path to the snp file.
#' @param ind_file Path to the ind file.
#' @param badsnp_file SNP file with information about ignored sites.
#' @export
f4_ratio <- function(X, A, B, C, O,
                     eigenstrat_prefix=NULL, geno=NULL, snp=NULL, ind=NULL, badsnp=NULL,
                     dir_name=NULL) {
    # get the path to the population, parameter and log files
    prefix <- paste0("f4_ratio_", A, "_", B, "_", C, "_", O)
    files <- get_files(dir_name, prefix)

    create_qpF4ratio_pops(X=samples$name, A=A, B=B, C=C, O=O, file=files[["pop_file"]])
    create_param_file(files[["par_file"]], files[["pop_file"]], eigenstrat_prefix,
                      geno, snp, ind, badsnp)

    run_cmd("qpF4ratio", param_file=files[["param_file"]], log_file=files[["log_file"]])
    
    read_qpF4ratio(paste0(prefix, ".log")) %>% mutate(setup=prefix)
}


# Reading output log files ====================================================


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



# Generating config files ==================================================


#' Generate a file with population "constellations" for a qpF4ratio run.
#'
#' @description For details about the format of the population list
#'     file see documentation of the ADMIXTOOLS package.
#' 
#' @param X, A, B, C, O Population names, using the terminology of
#'     Patterson et al., 2012
#' @param file Path to the poplist file that will be generated.
#' @export
create_qpF4ratio_pop_file <- function(X, A, B, C, O, file) {
    lines <- sprintf("%s %s : %s %s :: %s %s : %s %s", A, O, X, C, A, O, B, C)
    writeLines(lines, file)
}


#' Generate a file with population "constellations" for a qpF4ratio run.
#'
#' @description For details about the format of the population list
#'     file see documentation of the ADMIXTOOLS package.
#' 
#' @param X, A, B, C, O Population names, using the terminology of
#'     Patterson et al., 2012. It is possible to specify vectors of
#'     population identifiers as W, X, etc. That way, all possible
#'     combinations of specified W, X, Y and Z quadruples will be
#'     generated.
#' @param poplist List of populations. It will be saved as a file, one
#'     population per line.
#' @param file Path to the poplist file that will be generated.
#' @export
create_Dstats_pop_file <- function(W=NULL, X=NULL, Y=NULL, Z=NULL, poplist=NULL, file) {    
    if (all(!is.null(c(W, X, Y, Z)))) { # individual populations were specified
        lines <- c()
        for (w in W) for (x in X) for (y in Y) for (z in Z) {
            lines <- c(lines, sprintf("%s %s %s %s", w, x, y, z))
        }
    } else if (!is.null(poplist)) { # a simple list of populations was specified
        lines <- poplist
    } else { # missing population information
        stop("Population identifiers must be specified")
    }
    
    writeLines(lines, file)
}


#' Generate a parameter file.
#' 
#' @param param_file Path to the parameter file generated by this
#'     function.
#' @param pop_file Populations file (qpDstat quadruples, Dstats
#'     population list, qpF4ratio populations, etc.).
#' @param eigenstrat_prefix Prefix of the geno/snp/ind files (can
#'     include the path). If specified, geno_file/snp_file/ind_file
#'     will be ignored.
#' @param geno_file Path to the genotype file.
#' @param snp_file Path to the snp file.
#' @param ind_file Path to the ind file.
#' @param badsnp_file SNP file with information about ignored sites.
#' @param f4mode Run Dstats with an "f4mode: YES" option, estimating
#'     just the f4?
#' @export
create_param_file <- function(param_file, pop_file,
                              eigenstrat_prefix=NULL,
                              geno_file=NULL, snp_file=NULL, ind_file= NULL,
                              badsnp_file=NULL,
                              f4mode=FALSE) {
    if (!is.null(eigenstrat_prefix)) {
        geno_file <- paste0(eigenstrat_prefix, ".geno")
        snp_file <- paste0(eigenstrat_prefix, ".snp")
        ind_file <- paste0(eigenstrat_prefix, ".ind")
    }

    if (any(is.null(c(geno_file, snp_file, ind_file)))) {
        stop("Either the 'eigenstrat_prefix' or the paths to geno/snp/ind files must be specified")
    }

    writeLines(sprintf("genotypename: %s\nsnpname: %s\nindivname: %s\npopfilename: %s",
                       geno_file, snp_file, ind_file, pop_file),
               con=param_file)

    if (!is.null(badsnp_file)) {
        write(sprintf("badsnpname: %s", badsnp_file), file=param_file, append=TRUE)
    }

    if (f4mode) {
        write("f4mode: YES", file=param_file, append=TRUE)
    }
}



# Running ADMIXTOOLS commands =================================================


#' Run a given ADMIXTOOLS command.
#'
#' @param tool Which ADMIXTOOLS command to run. At the moment, only
#'     qpDstat and qpF4ratio are supported.
#' @param param_file Path to the parameter file.
#' @param log_file Path to the output file. If NULL, output will be
#'     printed to stdout.
#' @export
run_cmd <- function(cmd, param_file, log_file) {
    if (!cmd %in% c("qpDstat", "qpF4ratio")) {
        stop("ADMIXTOOLS command '", cmd, "' is not supported or does not exist")
    }

    output <- paste("> ", log_file)

    system(paste(cmd, "-p", param_file, output))
}



# EIGENSTRAT manipulation utilities ============================================


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
    read_table(file, col_names=c("id", "sex", "label"))
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



# Utility functions for filtering  ============================================


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

