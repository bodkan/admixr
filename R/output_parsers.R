#' Read output log file from a qpF4ratio run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with parsed results.
#'
#' @export
#'
#' @importFrom magrittr "%>%"
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
#' @return Tibble object with parsed results.
#'
#' @export
#'
#' @importFrom magrittr "%>%"
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
      raw_cols[, !purrr::map_lgl(raw_cols,
                                 ~ any(stringr::str_detect(., "best")))] %>%
      {
        setNames(., c("W", "X", "Y", "Z", result_col,
                      if (ncol(.) > 9) { "stderr" } else{ NULL },
                     "Zscore", "BABA", "ABBA", "n_snps"))
      }

    res_df
}


#' Read output log file from a qp3Pop run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with parsed results.
#'
#' @export
#'
#' @importFrom magrittr "%>%"
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
        setNames(c("A", "B", "C", "f3", "stderr", "Zscore", "n_snps"))

    res_df
}
