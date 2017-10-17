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
    setup <- paste0("qpDstat")
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qpDstat_pop_file(W, X, Y, Z, file=files[["pop_file"]])
    create_par_file(files[["par_file"]], files[["pop_file"]],
                    prefix, geno, snp, ind, badsnp, f4mode)

    if (f4mode) {
        write("f4mode: YES", file=files[["par_file"]], append=TRUE)
    }

    run_cmd("qpDstat", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qpDstat(files[["log_file"]])
}



#' Calculate the 3-population statistic and return the results as a data.frame.
#'
#' @param A, B, C Population names, using the terminology of Patterson
#'     et al., 2012
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix'
#'     argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @export
qp3Pop <- function(A, B, C,
                   prefix=NULL, geno=NULL, snp=NULL, ind=NULL, badsnp=NULL,
                   dir_name=NULL, inbreed=FALSE) {
    check_presence(c(A, B, C), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qp3Pop")
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qp3Pop_pop_file(A, B, C, file=files[["pop_file"]])
    create_par_file(files[["par_file"]], files[["pop_file"]],
                    prefix, geno, snp, ind, badsnp, f4mode=FALSE, inbreed)

    if (inbreed) {
        write("inbreed: YES", file=files[["par_file"]], append=TRUE)
    }

    run_cmd("qp3Pop", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qp3Pop(files[["log_file"]])
}



#' Run the qpAdm analysis.
#'
#' @param L, R Sets of left (U) and right (R) populations using the
#'     terminology of Haak et al., 2012 (Supplementary Information 10
#'     on page 128).
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix'
#'     argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @export
qpAdm <- function(L, R,
                  prefix=NULL, geno=NULL, snp=NULL, ind=NULL, badsnp=NULL,
                  dir_name=NULL) {
    check_presence(c(L, R), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qpAdm")
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qpAdm_pop_files(L, R, config_prefix)
    create_par_file(files[["par_file"]], files[["pop_file"]],
                    prefix, geno, snp, ind, badsnp)

    run_cmd("qpAdm", par_file=files[["par_file"]], log_file=files[["log_file"]])

#    read_qp3Pop(files[["log_file"]])
}
