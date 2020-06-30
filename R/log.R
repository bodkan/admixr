#' Print the full log output of an admixr wrapper to the console.
#'
#' @param x Output from one of the admixr wrappers (d, f4, qpAdm, ...)
#' @param target A specific log to examine (relevant for multiple target qpAdm runs)
#' @param save Save the log output to a disk?
#' @param prefix Prefix of the output log file(s) (name of the admixr command by default)
#' @param dir In which directory to save the log file(s)?
#' @param suffix Suffix of the output log file(s) (".txt" by default)
#'
#' @examples
#' \dontrun{# download an example genomic data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # define a set of populations to analyze and calculate a D statistic
#' pops <- c("French", "Sardinian", "Han", "Papuan", "Khomani_San", "Mbuti", "Dinka")
#' result_d <- d(
#'     W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp",
#'     data = snps
#' )
#'
#' # examine the full log output associated with the returned object
#' loginfo(result_d)
#' }
#'
#' @export
loginfo <- function(x, target = NA, save = FALSE, prefix = NA, dir = ".", suffix = ".txt") {
  if (!inherits(x, "admixr_result")) {
    stop("Object does not contain a result of an admixr function", call. = FALSE)
  }

  # extract the hidden admixr attributes for easier handling
  cmd <- attr(x, "command")
  log_output <- attr(x, "log_output")
  targets <- names(log_output)

  if (!is.na(target) && !cmd %in% c("qpAdm", "qpAdm_rotation"))
    stop(glue::glue("Specifying target does not make sense for examining the log output of {cmd}"),
         call. = FALSE)

  if (!is.na(target) && !target %in% targets)
    stop(glue::glue("Target/model '{target}' is not present in the output (choices are: {paste(targets, collapse = ', ')})"),
         call. = FALSE)

  ## qpAdm's log output are stored as a list of character vectors but everything
  ## else is simply a character vector - we convert everything to a list to
  ## iterate over log outputs below
  if (!is.list(log_output))
    log_output <- list(log_output)

  for (i in seq_along(log_output)) {
      ## write only single target/model log if requested by the user
      if (!is.na(target) && !is.null(targets) && target != targets[i]) next

    if (save) {
      if (is.na(prefix)) prefix <- cmd
  
      # generate output file name
      if (length(log_output) > 1)
        output <- glue::glue("{prefix}_{targets[i]}{suffix}")
      else
        output <- glue::glue("{prefix}{suffix}")

      writeLines(log_output[[i]], con = file.path(dir, output))
    } else {
      if (cmd == "qpAdm") {
        title <- glue::glue("qpAdm for target '{targets[i]}'")
      } else if (cmd == "qpAdm_rotation") {
        title <- glue::glue("qpAdm rotation for model '{targets[i]}'")
      } else {
        title <- cmd
      }

      if (i > 1 && is.na(target)) cat("\n\n")
      cat(paste0("Full output log of ", title, ":\n"))
      cat("===================================================\n\n")
      cat(paste(log_output[[i]], collapse = "\n"))
      cat("\n")
    }
  }
}


#' Print out the admixr result object (dataframe or a list) without showing
#' the hidden attributes.
#'
#' @param x admixr output object (dataframe or a list produced by qpAdm/qpWave)
#' @param ... Additional arguments passed to print.
#'
#' @export
print.admixr_result <- function(x, ...) {
    if (attr(x, "command") %in% c("qpAdm", "qpAdm_rotation") && length(x) == 3) {
    print.default(list(
      proportions = x$proportions,
      ranks = x$ranks,
      subsets = x$subsets
    ), ...)
  } else {
    print(tibble::as_tibble(x))
  }
}

