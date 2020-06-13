#' Print full output log associated with an admixr output.
#'
#' @param x Output from one of the admixr wrappers
#' @param target Specify a specific log to examine (relevant multiple target qpAdm run)
#' @param save Save the log output(s) to a file?
#' @param prefix Path prefix to the saved log output (name of the admixr command by default)
#'
#' @export
printlog <- function(x, target = NA, save = FALSE, prefix = NA) {
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
        output <- glue::glue("{prefix}_{targets[i]}.txt")
      else
        output <- glue::glue("{prefix}.txt")

      writeLines(log_output[[i]], con = output)
    } else {
      if (cmd == "qpAdm") {
        title <- glue::glue("qpAdm for target '{targets[i]}'")
      } else {
        title <- cmd
      }
      
      if (i > 1) cat("\n\n")
      cat(paste0("Full output log of ", title, ":\n"))
      cat("==================================================\n\n")
      cat(paste(log_output[[i]], collapse = "\n"))
    }
  }
}


print.admixr_result <- function(x, ...) {
  if (attr(x, "command") == "qpAdm") {
    print.default(list(
      proportions = x$proportions,
      ranks = x$ranks,
      subsets = x$subsets
    ), ...)
  } else {
    print(tibble::as_tibble(x, ...))
  }
}

