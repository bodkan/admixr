#' Print the full log output associated with an admixr output object.
#'
#' @param x Output from one of the admixr wrappers (d, f4, qpAdm, ...)
#' @param target A specific log to examine (relevant for multiple target qpAdm runs)
#' @param save Save the log output to a disk?
#' @param prefix Prefix of the output log file(s) (name of the admixr command by default)
#' @param dir In which directory to save the log file(s)?
#' @param suffix Suffix of the output log file(s) (".txt" by default)
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

  if (!is.na(target) && cmd != "qpAdm")
    stop(glue::glue("Specifying target does not make sense for examining the log output of {cmd}"),
         call. = FALSE)

  if (!is.na(target) && !target %in% targets)
    stop(glue::glue("Target '{target}' is not present in the output (choices are: {paste(targets, collapse = ', ')})"),
         call. = FALSE)

  # qpAdm's log output are stored as a list of character vectors but everything
  # else is simply a character vector - we convert everything to a list to
  # iterate over log outputs below
  if (!is.list(log_output))
    log_output <- list(log_output)

  for (i in seq_along(log_output)) {
    # write only single target log if requested by the user
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
      } else {
        title <- cmd
      }
      
      if (i > 1 && is.na(target)) cat("\n\n")
      cat(paste0("Full output log of ", title, ":\n"))
      cat("==================================================\n\n")
      cat(paste(log_output[[i]], collapse = "\n"))
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
  if (attr(x, "command") == "qpAdm") {
    print.default(list(
      proportions = x$proportions,
      ranks = x$ranks,
      subsets = x$subsets
    ), ...)
  } else {
    print(tibble::as_tibble(x))
  }
}

