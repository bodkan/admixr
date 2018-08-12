#' Read output file of one of the ADMIXTOOLS' programs.
#' 
#' This function servers as a dispatcher delegating the parsing to
#' read_qpF4ratio, read_qpDstat, etc.
#'
#' @param file Name of the output log file.
#'
#' @return data.frame with parsed results.
#'
#' @export
read_output <- function(file) {
  # extract ADMIXTOOLS command name from the output file
  cmd <- readLines(file) %>%
    stringr::str_match("^## (\\w+) version:") %>%
    .[complete.cases(.)] %>% .[2]
  
  parsers <- list(
    qp3Pop = read_qp3Pop,
    qpDstat = read_qpDstat,
    qpF4ratio = read_qpF4ratio,
    qpAdm = read_qpAdm
  )
  
  # it feels a little dumb, re-reading the whole output file a 2nd time,
  # there must be a cleaner way to do this
  as.data.frame(parsers[[cmd]](file))
}


#' Read output log file from a qpF4ratio run.
#'
#' @param file Name of the output log file.
#'
#' @return data.frame with parsed results.
read_qpF4ratio <- function(file) {
    log_lines <- readLines(file) %>% .[!stringr::str_detect(., "warning")]

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
        setNames(c("A", "O", "X", "C", "A", "O", "B", "C", "alpha", "stderr", "Zscore")) %>%
        .[c("A", "B", "X", "C", "O", "alpha", "stderr", "Zscore")]

    res_df
}


#' Read output log file from a qpDstat run.
#'
#' @param file Name of the output log file.
#'
#' @return data.frame with parsed results.
read_qpDstat <- function(file) {
    log_lines <- readLines(file) %>%
      .[!stringr::str_detect(., "warning")] %>%
      .[!stringr::str_detect(., "nodata")]

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
      raw_cols[, !sapply(raw_cols,
                         function(col) any(stringr::str_detect(col, "best")))] %>%
      {
        setNames(., c("W", "X", "Y", "Z", result_col,
                      if (ncol(.) > 9) { "stderr" } else{ NULL },
                     "Zscore", "BABA", "ABBA", "nsnps"))
      }

    res_df
}


#' Read output log file from a qp3Pop run.
#'
#' @param file Name of the output log file.
#'
#' @return data.frame with parsed results.
read_qp3Pop <- function(file) {
    log_lines <- readLines(file)

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[stringr::str_detect(log_lines, "result:")] %>%
        stringr::str_replace("result: ", "") %>%
        stringr::str_replace_all(" +", " ") %>%
        stringr::str_replace_all("^ | $", "")

    res_df <- res_lines %>%
        paste0("\n", collapse = "\n") %>%
        readr::read_delim(delim = " ", col_names = FALSE) %>%
        setNames(c("A", "B", "C", "f3", "stderr", "Zscore", "nsnps"))

    res_df
}


#' Read output log file from a qp3Pop run.
#'
#' @param file Name of the output log file.
#'
#' @return data.frame object with parsed results.
read_qpAdm <- function(file) {
  log_lines <- readLines(file)

  # parse the lines of the results section and extract the names of
  # tested populations/individuals, estimated admixture proportions
  # alpha, std. errors and Z-score
  stats <- stringr::str_subset(log_lines, "(Jackknife mean|std. errors):") %>%
    stringr::str_replace("(Jackknife mean|std. errors): +", "") %>%
    stringr::str_replace_all(" +", " ") %>%
    stringr::str_replace_all("^ | $", "") %>%
    stringr::str_split(" ") %>%
    lapply(as.numeric) %>%
    setNames(c("proportion", "stderr"))
  leftpops <- stringr::str_locate(log_lines, "(left|right) pops:") %>%
    .[, 1] %>%  { !is.na(.) } %>% which

  target_pop <- log_lines[leftpops[1] + 1]
  source_pops <- log_lines[(leftpops[1] + 2) : (leftpops[2] - 2)]

  snp_count <- stringr::str_subset(log_lines, paste0("coverage: +", target_pop)) %>%
    stringr::str_replace(paste0("coverage: +", target_pop, " +"), "") %>%
    as.numeric

  # wide format
  rbind(c(target_pop, snp_count, stats$proportion, stats$stderr)) %>%
    tibble::as_tibble() %>%
    setNames(c("target", "nsnps", source_pops, paste0("stderr_", source_pops))) %>%
    dplyr::mutate_at(dplyr::vars(-target), as.numeric)

  # # long format
  # tibble::tibble(
  #   target_pop,
  #   snp_count,
  #   source_pops,
  #   proportion = stats$proportion,
  #   stderr = stats$stderr
  # )
}
