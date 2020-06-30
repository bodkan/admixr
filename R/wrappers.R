#' Calculate the D, f4, f4-ratio, or f3 statistic.
#' @param W,X,Y,Z,A,B,C,O Population names according to the nomenclature used in
#'     Patterson et al., 2012.
#'
#' @inheritParams qpAdm
#'
#' @return Data frame object with calculated statistics
#'
#' @examples
#' \dontrun{# download an example genomic data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # define a set of populations to analyze
#' pops <- c("French", "Sardinian", "Han", "Papuan", "Dinka")
#'
#' result_f4ratio <- f4ratio(
#'     X = pops, A = "Altai", B = "Vindija", C = "Yoruba", O = "Chimp",
#'     data = snps
#' )
#'
#' result_d <- d(
#'     W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp",
#'     data = snps
#' )
#'
#' result_f4 <- f4(
#'     W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp",
#'     data = snps
#' )
#'
#' result_f3 <- f3(
#'     A = pops, B = "Mbuti", C = "Khomani_San",
#'     data = snps
#' )
#' }
#'
#' @export
f4ratio <- function(data, X, A, B, C, O, outdir = NULL, params = NULL) {
  check_presence(c(X, A, B, C, O), data)

  # get the path to the population, parameter and log files
  config_prefix <- paste0("qpF4ratio__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
  files <- get_files(outdir, config_prefix)

  create_qpF4ratio_pop_file(X = X, A = A, B = B, C = C, O = O, file = files[["pop_file"]])
  create_par_file(files, data, params)

  run_cmd("qpF4ratio", par_file = files[["par_file"]], log_file = files[["log_file"]])

  read_output(files[["log_file"]])
}



#' @rdname f4ratio
#'
#' @param f4mode Calculate the f4 statistic instead of the D statistic.
#' @param quartets List of character vectors (quartets of population/sample labels)
#'
#' @export
d <- function(data, W, X, Y, Z, quartets = NULL, outdir = NULL, f4mode = FALSE, params = NULL) {
  if (is.null(quartets)) {
    check_presence(c(W, X, Y, Z), data)
  } else {
    check_presence(unique(unlist(quartets)), data)
  }

  # get the path to the population, parameter and log files
  config_prefix <- paste0("qpDstat__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
  files <- get_files(outdir, config_prefix)

  if (is.null(quartets)) {
    create_qpDstat_pop_file(W, X, Y, Z, file = files[["pop_file"]])
  } else {
    create_qpDstat_pop_file_quartets(quartets, file = files[["pop_file"]])
  }
  create_par_file(files, data, params)

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
f4 <- function(data, W, X, Y, Z, quartets = NULL, outdir = NULL, params = NULL) {
  d(data, W, X, Y, Z, quartets, outdir, f4mode = TRUE)
}



#' @rdname f4ratio
#'
#' @param inbreed See README.3PopTest in ADMIXTOOLS for an explanation.
#'
#' @export
f3 <- function(data, A, B, C, outdir = NULL, inbreed = FALSE, params = NULL) {
  check_presence(c(A, B, C), data)

  # get the path to the population, parameter and log files
  config_prefix <- paste0("qp3Pop__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
  files <- get_files(outdir, config_prefix)

  create_qp3Pop_pop_file(A, B, C, file = files[["pop_file"]])
  create_par_file(files, data, params)

  if (inbreed) {
    write("inbreed: YES", file = files[["par_file"]], append = TRUE)
  }

  run_cmd("qp3Pop", par_file = files[["par_file"]], log_file = files[["log_file"]])

  read_output(files[["log_file"]])
}



#' Calculate ancestry proportions in a set of target populations.
#'
#' @param target Vector of target populations (evaluated one at a time).
#' @param sources Source populations related to true ancestors.
#' @param outgroups Outgroup populations.
#'
#' @param data EIGENSTRAT data object.
#' @param outdir Where to put all generated files (temporary directory by default).
#' @param params Named list of parameters and their values.
#'
#' @return List of three components:
#'     1. estimated ancestry proportions
#'     2. ranks statistics
#'     3. analysis of patterns (all possible subsets of ancestry sources).
#'
#' @examples
#' \dontrun{# download example data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # estimate the proportion of Neandertal ancestry in a French
#' # individual and other associated qpAdm statistics (see detailed
#' # description in the tutorial vignette)
#' result <- qpAdm(
#'     target = "French",
#'     sources = c("Vindija", "Yoruba"),
#'     outgroups = c("Chimp", "Denisova", "Altai"),
#'     data = snps
#' )
#' }
#'
#' @export
qpAdm <- function(data, target, sources, outgroups, outdir = NULL,
                  params = list(allsnps = "YES", summary = "YES", details = "YES")) {
    if (length(outgroups) < length(sources) + 1) {
        stop("The number of outgroup samples has to be larger or equal than the number of sources + 1",
             call. = FALSE)
    }

  check_presence(c(target, sources, outgroups), data)

  results <- lapply(target, function(X) {
    # get the path to the population, parameter and log files
    config_prefix <- paste0("qpAdm__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
    files <- get_files(outdir, config_prefix)

    files[["popleft"]] <-  stringr::str_replace(files[["pop_file"]], "$", "left")
    files[["popright"]] <-  stringr::str_replace(files[["pop_file"]], "$", "right")
    files[["pop_file"]] <- NULL

    create_leftright_pop_files(c(X, sources), outgroups, files)
    create_par_file(files, data, params)

    run_cmd("qpAdm", par_file = files[["par_file"]], log_file = files[["log_file"]])

    read_output(files[["log_file"]])
  })

  # process the complex list of lists of dataframes into a more readable form
  # by concatenating all internal dataframes and returning a simple list
  # of three dataframes
  proportions <- dplyr::bind_rows(lapply(results, function(x) x$proportions))

  ranks <- lapply(seq_along(target), function(i) { results[[i]]$ranks %>% dplyr::mutate(target = target[i]) }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(target, dplyr::everything())

  res <- list(
    proportions = proportions,
    ranks = ranks
  )

  if (!is.null(results[[1]]$subsets)) {
    subsets <- lapply(seq_along(target), function(i) {
        results[[i]]$subsets %>% dplyr::mutate(target = target[i])
      }) %>%
      dplyr::bind_rows() %>%
      dplyr::select(target, dplyr::everything())
    res$subsets <- subsets
  } else {
    res$subsets <- NULL
  }
  
  attr(res, "command") <- "qpAdm"
  attr(res, "log_output") <- lapply(results, function(i) attr(i, "log_output"))
  # name the elements of the list of log outputs based on target population
  names(attr(res, "log_output")) <- target
  class(res) <- c("admixr_result", class(res))

  res
}



#' Find the most likely number of ancestry waves using the qpWave method.
#'
#' Given a set of 'left' populations, estimate the lowest number of necessary
#' admixture sources related to the set of 'right' populations.
#'
#' It has been shown (Reich, Nature 2012 - Reconstructing Native American
#' population history) that if the 'left' populations are mixtures of N
#' different sources related to the set of 'right' populations, the rank of the
#' matrix of the form \eqn{f_4(left_i, left_j; right_k, right_l)} will have a
#' rank N - 1. This function uses the ADMIXTOOLS command qpWave to find the
#' lowest possible rank of this matrix that is consistent with the data.
#'
#' @param left,right Character vectors of populations labels.
#' @param maxrank Maximum rank to test for.
#' @param details Return the A, B matrices used in rank calculations?
#' @inheritParams qpAdm
#'
#' @return Table of rank test results.
#'
#' @examples
#' \dontrun{# download example data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # run the qpWave wrapper (detailed description in the tutorial vignette)
#' result <- qpWave(
#'      left = c("French", "Sardinian", "Han"),
#'      right = c("Altai", "Yoruba", "Mbuti"),
#'      data = snps
#' )
#' }
#'
#' @export
qpWave <- function(data, left, right, maxrank = NULL, details = FALSE, outdir = NULL, params = NULL) {
  check_presence(c(left, right), data)
  if (length(intersect(left, right))) {
    stop("Duplicated populations in both left and right population sets not allowed: ",
         paste(intersect(left, right), collapse = " "),
         call. = FALSE)
  }

  # get the path to the population, parameter and log files
  setup <- paste0("qpWave")
  config_prefix <- paste0(setup, "__", as.integer(stats::runif(1, 0, .Machine$integer.max)))
  files <- get_files(outdir, config_prefix)

  files[["popleft"]] <-  stringr::str_replace(files[["pop_file"]], "$", "left")
  files[["popright"]] <-  stringr::str_replace(files[["pop_file"]], "$", "right")
  files[["pop_file"]] <- NULL

  create_leftright_pop_files(left, right, files)
  create_par_file(files, data, params)

  if (!is.null(maxrank)) {
    write(sprintf("maxrank: %d", maxrank), file = files[["par_file"]], append = TRUE)
  }

  run_cmd("qpWave", par_file = files[["par_file"]], log_file = files[["log_file"]])

  read_output(files[["log_file"]], details)
}
