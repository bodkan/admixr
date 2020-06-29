#' Read an output file from one of the ADMIXTOOLS programs.
#'
#'@param file A path to an output file.
#'@param ... See the 'details' argument of qpWave.
#'
#'@return A tibble with the results.
#'
#'@export
read_output <- function(file, ...) {
  log_lines <- readLines(file) %>% .[!stringr::str_detect(., "nodata")]

  # Sometimes ADMIXTOOLS crashes for various reasons, leaving a truncated
  # log file - this leads to cryptic nonsensical error messages on our part,
  # which only confuse the user. Let's make sure that the log file is
  # complete, informing the user to inspect it in case we see a problem.
  if (!any(stringr::str_detect(log_lines, "end of"))) {
    message("BEGINNING OF OUTPUT FILE\n",
            "==================================================\n\n",
            paste(log_lines, collapse = "\n"),
            "\n\n==================================================\n",
            "END OF OUTPUT FILE\n\n")
    stop("The output file we got from ADMIXTOOLS is truncated.
       Please examine the full output above and check for errors.",
         call. = FALSE)
  }

  # extract ADMIXTOOLS command name from the output file
  cmd <- log_lines %>%
    stringr::str_match("^## (\\w+) version:") %>%
    .[stats::complete.cases(.)] %>% .[2]

  parsers <- list(
    qp3Pop = read_qp3Pop,
    qpDstat = read_qpDstat,
    qpF4ratio = read_qpF4ratio,
    qpWave = read_qpWave,
    qpAdm = read_qpAdm
  )

  # it feels a little dumb, re-reading the whole output file a 2nd time,
  # there must be a cleaner way to do this
  result <- parsers[[cmd]](log_lines, ...)

  # annotate the object with metadata and the full log output
  attr(result, "command") <- cmd
  attr(result, "log_output") <- log_lines
  class(result) <- c("admixr_result", class(result))

  result
}



# Read output log file from a qpF4ratio run.
read_qpF4ratio <- function(log_lines) {
  # extract the number of analyzed test populations/individuals
  # (corresponding to the number of rows of the results table)
  n_pops <- log_lines[which(stringr::str_detect(log_lines, "^nplist:"))] %>%
    stringr::str_extract("[0-9]+$") %>%
    as.integer

  # parse the lines of the results section and extract the names of
  # tested populations/individuals, estimated admixture proportions
  # alpha, std. errors and Z-score
  res_lines <- log_lines[(length(log_lines) - n_pops) : (length(log_lines) - 1)] %>%
    stringr::str_replace("result: ", "") %>%
    stringr::str_replace_all(":", "") %>%
    stringr::str_replace_all(" +", " ") %>%
    stringr::str_replace("^ ", "")

  res_df <- res_lines %>%
    paste0("\n", collapse = "\n") %>%
    readr::read_delim(delim = " ", col_names = FALSE) %>%
    stats::setNames(c("A", "O", "X", "C", "A", "O", "B", "C", "alpha", "stderr", "Zscore")) %>%
    .[c("A", "B", "X", "C", "O", "alpha", "stderr", "Zscore")]

  res_df
}



# Read output log file from a qpDstat run.
read_qpDstat <- function(log_lines) {
  # check for duplicated populations:
  #   - if there is no "result" row, issue a warning
  #   - otherwise ignore the "repeated populations" warning
  if (!any(stringr::str_detect(log_lines, "result: "))) {
    warning("ADMIXTOOLS does not allow duplicated population names when calculating a single f4 or D statistic.\nReturning an empty dataframe...",
            call. = FALSE)
    return(tibble::tibble())
  } else {
    log_lines <- log_lines %>% .[!stringr::str_detect(., "warning")]
  }

  # extract the number of analyzed population quadruples
  n_quads <- length(log_lines) - (which(stringr::str_detect(log_lines, "^nrows, ncols:"))) - 1

  # parse the lines of the results section and extract the names of
  # tested populations/individuals, estimated admixture proportions
  # alpha, std. errors and Z-score
  res_lines <- log_lines[(length(log_lines) - n_quads) : (length(log_lines) - 1)] %>%
    stringr::str_replace("result: ", "") %>%
    stringr::str_replace_all(" +", " ") %>%
    stringr::str_replace_all("^ | $", "")

  result_col <- ifelse(any(stringr::str_detect(log_lines, "f4mode: YES")), "f4", "D")

  raw_cols <- res_lines %>%
    paste0("\n", collapse = "\n") %>%
    readr::read_delim(delim = " ", col_names = FALSE)

  # remove the weird "best" column first, then add an optional stderr column
  # (if it's present)
  res_df <-
    raw_cols[, !sapply(raw_cols, function(col) any(stringr::str_detect(col, "best")))] %>%
    {
      stats::setNames(., c("W", "X", "Y", "Z", result_col,
                           if (ncol(.) > 9) { "stderr" } else{ NULL },
                           "Zscore", "BABA", "ABBA", "nsnps"))
    }

  res_df
}



# Read output log file from a qp3Pop run.
read_qp3Pop <- function(log_lines) {
  # parse the lines of the results section and extract the names of
  # tested populations/individuals, estimated admixture proportions
  # alpha, std. errors and Z-score
  res_lines <- log_lines[stringr::str_detect(log_lines, "result:|no data")] %>%
    stringr::str_replace("result: |no data", "") %>%
    stringr::str_replace_all(" +", " ") %>%
    stringr::str_replace_all("^ | $", "")

  res_df <- res_lines %>%
    paste0("\n", collapse = "\n") %>%
    readr::read_delim(delim = " ", col_names = FALSE) %>%
    stats::setNames(c("A", "B", "C", "f3", "stderr", "Zscore", "nsnps"))

  res_df
}



parse_matrix <- function(lines) {
  m <- lines %>%
    stringr::str_replace_all("^ +| +$", "") %>%
    stringr::str_replace_all(" +", "\t") %>%
    paste0(collapse = "\n") %>%
    readr::read_tsv(col_names = FALSE) %>%
    t

  # rename columns and reset rownames
  colnames(m) <- unlist(m[1, ])
  m <- m[-1, , drop = FALSE]
  rownames(m) <- NULL
  class(m) <- "numeric"

  m
}



# Read output log file from a qp3Pop run.
read_qpWave <- function(log_lines, details = FALSE) {
  test_pos  <- which(stringr::str_detect(log_lines, "f4info:"))
  b_pos <- which(stringr::str_detect(log_lines, "B:"))
  a_pos <- which(stringr::str_detect(log_lines, "A:"))
  a_end <- c(test_pos[-c(1, 2)], which(stringr::str_detect(log_lines, "## end of run")))

  test_df <- log_lines[test_pos + 1] %>%
    stringr::str_replace_all("[a-z0-9]+: ", "") %>%
    stringr::str_replace_all(" +", "\t") %>%
    paste0(collapse = "\n") %>%
    readr::read_tsv(col_names = c("rank", "df", "chisq", "tail", "dfdiff",
                                  "chisqdiff", "taildiff"))

  if (details) {
    B_matrix <- lapply(seq_along(b_pos), function(i) {
      parse_matrix(log_lines[(b_pos[i] + 1) : (a_pos[i] - 1)])
    }) %>% stats::setNames(paste0(seq_along(.)))
  
    A_matrix <- lapply(seq_along(a_pos), function(i) {
      parse_matrix(log_lines[(a_pos[i] + 1) : (a_end[i] - 2)])
    }) %>% stats::setNames(paste0(seq_along(.)))
  
    matrices <- lapply(seq_along(B_matrix), function(rank) {
      list(A = A_matrix[[rank]], B = B_matrix[[rank]])
    })

    return(list(ranks = test_df, matrices = matrices))
  } else {
    return(test_df)
  }
    
}



# Read output log file from a qp3Pop run.
read_qpAdm <- function(log_lines) {
  # parse the admixture proportions and standard errors
  stats <- stringr::str_subset(log_lines, "(best coefficients|std. errors):") %>%
    stringr::str_replace("(best coefficients|std. errors): +", "") %>%
    stringr::str_replace_all(" +", " ") %>%
    stringr::str_replace_all("^ | $", "") %>%
    stringr::str_split(" ") %>%
    lapply(as.numeric) %>%
    stats::setNames(c("proportion", "stderr"))

  # parse the population names
  leftpops <- stringr::str_locate(log_lines, "(left|right) pops:") %>%
    .[, 1] %>%  { !is.na(.) } %>% which
  target <- log_lines[leftpops[1] + 1]
  source <- log_lines[(leftpops[1] + 2) : (leftpops[2] - 2)]

  # parse the SNP count
  nsnps_target <- stringr::str_subset(log_lines, paste0("coverage: +", target, " ")) %>%
    stringr::str_replace(paste0("coverage: +", target, " +"), "") %>%
    as.integer

  nsnps_used <- stringr::str_subset(log_lines, paste0("numsnps used: +")) %>%
    stringr::str_replace("numsnps used: +", "") %>%
    as.integer

  proportions <- rbind(c(target, stats$proportion, stats$stderr, nsnps_used, nsnps_target)) %>%
    tibble::as_tibble(.name_repair = "minimal") %>%
    stats::setNames(c("target", source, paste0("stderr_", source), "nsnps_used", "nsnps_target")) %>%
    dplyr::mutate_at(dplyr::vars(-target), as.numeric)

  # parse the population combination patterns into a data.frame
  # (but only if there are multiple sources!)
  if (any(stringr::str_detect(log_lines, "single source. terminating"))) {
    pat_df <- NULL
  } else {
    pat_start <- stringr::str_detect(log_lines, "fixed pat") %>% which
    pat_end <- stringr::str_detect(log_lines, "best pat") %>% which
    patterns <- log_lines[pat_start : (pat_end[1] - 1)] %>%
      stringr::str_replace(" fixed", "") %>%
      stringr::str_replace(" prob", "") %>%
      stringr::str_replace(" pattern", "") %>%
      stringr::str_replace_all(" +", " ") %>%
      stringr::str_replace_all("^ | $", "")
    pat_header <- c(strsplit(patterns[1], " ")[[1]], source) %>%
      stringr::str_replace("pat", "pattern")
    if (any(stringr::str_detect(patterns, "infeasible"))) {
      pat_header <- c(pat_header, "comment")
      patterns[-1] <- sapply(patterns[-1], USE.NAMES = FALSE, function(l)
        if (stringr::str_detect(l, "infeasible")) l else paste0(l, " -"))
    }
    pat_df <- patterns[-1] %>%
      paste0(collapse = "\n") %>%
      readr::read_delim(delim = " ", col_names = FALSE) %>%
      stats::setNames(pat_header)
  }

  # parse the rank test results
  ranks <- read_qpWave(log_lines)

  # add a p-value to the proportions dataframe
  proportions$pvalue <- pat_df[1, "tail", drop = TRUE]
  list(proportions = proportions, ranks = ranks, subsets = pat_df)
}
