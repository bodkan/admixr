#' Perform an f4-ratio calculation and return its result as a data.frame.
#'
#' @param X Population names, using the terminology of Patterson (2012).
#' @param A Population names, using the terminology of Patterson (2012).
#' @param B Population names, using the terminology of Patterson (2012).
#' @param C Population names, using the terminology of Patterson (2012).
#' @param O Population names, using the terminology of Patterson (2012).
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix' argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#'
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
    create_par_file(files, prefix, geno, snp, ind, badsnp)

    run_cmd("qpF4ratio", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qpF4ratio(files[["log_file"]]) %>% dplyr::mutate(setup=setup)
}


#' Calculate D statistics or F4 statistics (which is just the
#' numerator of a D statistic) and return results as a data.frame.
#'
#' @param W Population names, using the terminology of Patterson (2012).
#' @param X Population names, using the terminology of Patterson (2012).
#' @param Y Population names, using the terminology of Patterson (2012).
#' @param Z Population names, using the terminology of Patterson (2012).
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix' argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @param f4mode Calculate f4 statistics instead of D statistic.
#'
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
    create_par_file(files, prefix, geno, snp, ind, badsnp)

    if (f4mode) {
        write("f4mode: YES", file=files[["par_file"]], append=TRUE)
    }

    # automatically calculate standard errors
    write("printsd: YES", file = files[["par_file"]], append = TRUE)

    run_cmd("qpDstat", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qpDstat(files[["log_file"]])
}


#' Calculate a 3-population statistic and return results as a data.frame.
#'
#' @param A Population names, using the terminology of Patterson (2012).
#' @param B Population names, using the terminology of Patterson (2012).
#' @param C Population names, using the terminology of Patterson (2012).
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix'
#'     argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @param inbreed See README.3PopTest in ADMIXTOOLS.
#'
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
    create_par_file(files, prefix, geno, snp, ind, badsnp)

    if (inbreed) {
        write("inbreed: YES", file=files[["par_file"]], append=TRUE)
    }

    run_cmd("qp3Pop", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qp3Pop(files[["log_file"]])
}


#' Run a qpAdm analysis.
#'
#' @param L Sets of "left" populations using the terminology of Haak (2012,
#'     Supplementary Information 10 on page 128).
#' @param R Sets of "right" populations using the terminology of Haak (2012,
#'     Supplementary Information 10 on page 128).
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix'
#'     argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#'
#' @export
qpAdm <- function(L, R,
                  prefix=NULL, geno=NULL, snp=NULL, ind=NULL, badsnp=NULL,
                  dir_name=NULL) {
    check_presence(c(L, R), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qpAdm")
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    # ugh... fix this! also - rename this "files" thing, doesn't make any sense...
    files[["pop_file"]] <- NULL
    files[["popleft"]] <-  paste0(config_prefix, "__left")
    files[["popright"]] <-  paste0(config_prefix, "__right")

    create_qpAdm_pop_files(L, R, files)
    create_par_file(files, prefix, geno, snp, ind, badsnp)

    run_cmd("qpAdm", par_file=files[["par_file"]], log_file=files[["log_file"]])

#    read_qp3Pop(files[["log_file"]])
}