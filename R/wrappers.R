#' Perform an f4-ratio calculation and return its result as a data.frame.
#'
#' @param X,A,B,C,O Population names according to nomenclature used in
#'     Patterson et al., 2012.
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno,snp,ind Path to the geno/snp/ind file. Each overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#'
#' @export
f4_ratio <- function(X, A, B, C, O,
                     prefix = NULL, geno = NULL, snp = NULL, ind = NULL, badsnp = NULL,
                     dir_name = NULL) {
    check_presence(c(X, A, B, C, O), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qpF4ratio__", A, "_", B, "_", C, "_", O)
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qpF4ratio_pop_file(X = X, A = A, B = B, C = C, O = O, file = files[["pop_file"]])
    create_par_file(files, prefix, geno, snp, ind, badsnp)

    run_cmd("qpF4ratio", par_file = files[["par_file"]], log_file = files[["log_file"]])

    read_output(files[["log_file"]])
}


#' Calculate D statistics or F4 statistics (which is just the
#' numerator of a D statistic) and return results as a data.frame.
#'
#' @param W,X,Y,Z Population names according to nomenclature in Patterson
#'     et al., 2012.
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno,snp,ind Path to the geno/snp/ind file. Each overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @param f4mode Calculate f4 statistics instead of D statistic.
#'
#' @export
d <- function(W, X, Y, Z,
              prefix = NULL, geno = NULL, snp = NULL, ind = NULL, badsnp = NULL,
              dir_name = NULL, f4mode = FALSE) {
    check_presence(c(W, X, Y, Z), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qpDstat")
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qpDstat_pop_file(W, X, Y, Z, file = files[["pop_file"]])
    create_par_file(files, prefix, geno, snp, ind, badsnp)

    if (f4mode) {
        write("f4mode: YES", file = files[["par_file"]], append = TRUE)
    }

    # automatically calculate standard errors
    write("printsd: YES", file = files[["par_file"]], append = TRUE)

    run_cmd("qpDstat", par_file = files[["par_file"]], log_file = files[["log_file"]])

    read_output(files[["log_file"]])
}

#' Calculate a 3-population statistic and return results as a data.frame.
#'
#' @param A,B,C Population names.
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno,snp,ind Path to the geno/snp/ind file. Each overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @param inbreed See README.3PopTest in ADMIXTOOLS.
#'
#' @export
f3 <- function(A, B, C,
               prefix = NULL, geno = NULL, snp = NULL, ind = NULL, badsnp = NULL,
               dir_name = NULL, inbreed = FALSE) {
    check_presence(c(A, B, C), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qp3Pop")
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qp3Pop_pop_file(A, B, C, file = files[["pop_file"]])
    create_par_file(files, prefix, geno, snp, ind, badsnp)

    if (inbreed) {
        write("inbreed: YES", file = files[["par_file"]], append = TRUE)
    }

    run_cmd("qp3Pop", par_file = files[["par_file"]], log_file = files[["log_file"]])

    read_output(files[["log_file"]])
}

#' Calculate admixture proportions in a target population using qpAdm method.
#'
#' @param target Target population to estimate admixture proportions for.
#' @param source Source populations that are related to ancestors of target.
#' @param outgroup Outgroup populations,B,C Population names.
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno,snp,ind Path to the geno/snp/ind file. Each overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#'
#' @export
qpAdm <- function(target, source, outgroup,
                  prefix = NULL, geno = NULL, snp = NULL, ind = NULL, badsnp = NULL,
                  dir_name = NULL) {
  check_presence(c(target, source, outgroup), prefix, ind)
  
  # get the path to the population, parameter and log files
  setup <- paste0("qpAdm")
  config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
  files <- get_files(dir_name, config_prefix)

  files[["popleft"]] <-  stringr::str_replace(files[["pop_file"]], "$", "left")
  files[["popright"]] <-  stringr::str_replace(files[["pop_file"]], "$", "right")
  files[["pop_file"]] <- NULL

  create_qpAdm_pop_files(c(target, source), outgroup, files)
  create_par_file(files, prefix, geno, snp, ind, badsnp)
  
  run_cmd("qpAdm", par_file = files[["par_file"]], log_file = files[["log_file"]])
  
  read_output(files[["log_file"]])
}
