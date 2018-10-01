#' Calculate the D, f4, f4-ratio, or f3 statistic.
#' @param W,X,Y,Z,A,B,C,O Population names according to the nomenclature used in
#'     Patterson et al., 2012.
#'
#' @inheritParams qpAdm
#'
#' @export
f4ratio <- function(data, X, A, B, C, O, outdir = NULL) {
  check_presence(c(X, A, B, C, O), data)

  # get the path to the population, parameter and log files
  setup <- paste0("qpF4ratio__", A, "_", B, "_", C, "_", O)
  config_prefix <- paste0(setup, "__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
  files <- get_files(outdir, config_prefix)

  create_qpF4ratio_pop_file(X = X, A = A, B = B, C = C, O = O, file = files[["pop_file"]])
  create_par_file(files, data)

  run_cmd("qpF4ratio", par_file = files[["par_file"]], log_file = files[["log_file"]])

  read_output(files[["log_file"]])
}



#' @rdname f4ratio
#'
#' @param f4mode Calculate the f4 statistic instead of the D statistic.
#'
#' @export
d <- function(data, W, X, Y, Z, outdir = NULL, f4mode = FALSE) {
  check_presence(c(W, X, Y, Z), data)

  # get the path to the population, parameter and log files
  setup <- paste0("qpDstat")
  config_prefix <- paste0(setup, "__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
  files <- get_files(outdir, config_prefix)

  create_qpDstat_pop_file(W, X, Y, Z, file = files[["pop_file"]])
  create_par_file(files, data)

  if (f4mode) {
    write("f4mode: YES", file = files[["par_file"]], append = TRUE)
  }

  # automatically calculate standard errors
  write("printsd: YES", file = files[["par_file"]], append = TRUE)

  run_cmd("qpDstat", par_file = files[["par_file"]], log_file = files[["log_file"]])

  read_output(files[["log_file"]])
}



#' @rdname f4ratio
#'
#' @export
f4 <- function(data, W, X, Y, Z, outdir = NULL) {
  d(data, W, X, Y, Z, outdir, f4mode = TRUE)
}



#' @rdname f4ratio
#'
#' @param inbreed See README.3PopTest in ADMIXTOOLS for an explanation.
#'
#' @export
f3 <- function(data, A, B, C, outdir = NULL, inbreed = FALSE) {
  check_presence(c(A, B, C), data)

  # get the path to the population, parameter and log files
  setup <- paste0("qp3Pop")
  config_prefix <- paste0(setup, "__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
  files <- get_files(outdir, config_prefix)

  create_qp3Pop_pop_file(A, B, C, file = files[["pop_file"]])
  create_par_file(files, data)

  if (inbreed) {
    write("inbreed: YES", file = files[["par_file"]], append = TRUE)
  }

  run_cmd("qp3Pop", par_file = files[["par_file"]], log_file = files[["log_file"]])

  read_output(files[["log_file"]])
}



#' Calculate ancestry proportions in target populations using qpAdm.
#'
#' @param target Vector of target populations (evaluated one at a time).
#' @param references Reference source populations related to true ancestors.
#' @param outgroups Outgroup populations.
#'
#' @param data EIGENSTRAT data object.
#' @param outdir Where to put all generated files (temporary directory by default).
#'
#' @export
qpAdm <- function(data, target, references, outgroups, outdir = NULL) {
  check_presence(c(target, references, outgroups), data)

  dplyr::bind_rows(lapply(target, function(X) {
    # get the path to the population, parameter and log files
    setup <- paste0("qpAdm")
    config_prefix <- paste0(setup, "__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
    files <- get_files(outdir, config_prefix)

    files[["popleft"]] <-  stringr::str_replace(files[["pop_file"]], "$", "left")
    files[["popright"]] <-  stringr::str_replace(files[["pop_file"]], "$", "right")
    files[["pop_file"]] <- NULL

    create_leftright_pop_files(c(X, references), outgroups, files)
    create_par_file(files, data)

    run_cmd("qpAdm", par_file = files[["par_file"]], log_file = files[["log_file"]])

    read_output(files[["log_file"]])
  }))
}



#' Find the most likely number of ancestry waves.
#'
#' @param left,right Character vectors of populations labels.
#' @param maxrank Maximum rank to test for.
#' @param matrices Return the A, B matrices used in rank calculations?
#' @inheritParams qpAdm
#'
#' @export
qpWave <- function(data, left, right, maxrank = NULL, matrices = FALSE, outdir = NULL) {
  check_presence(c(left, right), data)

  # get the path to the population, parameter and log files
  setup <- paste0("qpWave")
  config_prefix <- paste0(setup, "__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
  files <- get_files(outdir, config_prefix)

  files[["popleft"]] <-  stringr::str_replace(files[["pop_file"]], "$", "left")
  files[["popright"]] <-  stringr::str_replace(files[["pop_file"]], "$", "right")
  files[["pop_file"]] <- NULL

  create_leftright_pop_files(left, right, files)
  create_par_file(files, data)

  if (!is.null(maxrank)) {
    write(sprintf("maxrank: %d", maxrank), file = files[["par_file"]], append = TRUE)
  }

  run_cmd("qpWave", par_file = files[["par_file"]], log_file = files[["log_file"]])

  read_output(files[["log_file"]], matrices)
}

