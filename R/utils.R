#' Count the number/proportion of present/missing sites in each sample
#'
#' @param data EIGENSTRAT data object.
#' @param prop Calculate the proportion instead of counts?
#' @param missing Count present SNPs or missing SNPs?
#'
#' @return A data.frame object with SNP counts/proportions.
#'
#' @examples
#' \dontrun{snps <- eigenstrat(download_data(dirname = tempdir()))
#' 
#' present_count <- count_snps(snps)
#' missing_count <- count_snps(snps, missing = TRUE)
#'
#' present_proportion <- count_snps(snps, prop = TRUE)
#' missing_proportion <- count_snps(snps, missing = TRUE, prop = TRUE)
#' }
#'
#' @export
#' @import rlang
count_snps <- function(data, missing = FALSE, prop = FALSE) {
  summary_fun <- ifelse(prop, mean, sum)
  if (missing) {
    fun <- function(x) as.integer(is.na(x))
    col <- "missing"
  } else {
    fun <- function(x) as.integer(!is.na(x))
    col <- "present"
  }
  geno <- read_geno(data)
  if (!is.null(data$exclude)) {
    snps <- read_snp(data)
    excluded_snps <- read_snp(data, exclude = TRUE)
    chroms <- unique(snps$chrom)
    sites_to_remove <- lapply(chroms, function(chrom) {
      chrom_snps <- snps[snps$chrom == chrom, ]
      chrom_excluded_snps <- excluded_snps[excluded_snps$chrom == chrom, ]
      chrom_snps$pos %in% chrom_excluded_snps$pos
    }) %>% Reduce(c, .)
    geno <- geno[!sites_to_remove, ]
  }
  result <- read_ind(data)
  result[[col]] <- as.vector(t(dplyr::summarise_all(geno, list(~ summary_fun(fun(.))))))
  result
}



# Run a specified ADMIXTOOLS command.
run_cmd <- function(cmd, par_file, log_file, directory = NULL) {
  # if needed, change into the required directory first
  cd_cmd <- if (!is.null(directory)) paste("cd", directory, ";") else ""
  system(paste(cd_cmd, cmd, "-p", par_file, ">", log_file))
}



# Create either specified or a temporary directory.
get_dir <- function(dir_name = NULL) {
  if (!is.null(dir_name)) {
    dir.create(dir_name, showWarnings = FALSE)
  } else {
    dir_name <- tempdir()
  }

  path.expand(dir_name)
}



# Generate paths to the population file, parameter file and log file
# based on a specified directory.
get_files <- function(dir_name, prefix) {
  directory <- get_dir(dir_name)
  list(
    pop_file = file.path(directory, paste0(prefix, ".pop")),
    par_file = file.path(directory, paste0(prefix, ".par")),
    log_file = file.path(directory, paste0(prefix, ".log"))
  )
}



# Check that the provided object is of the required type
check_type <- function(x, type) {
  if (!inherits(x, type)) {
    stop(glue::glue("Input variable '{deparse(substitute(x))}' is not of the type {type}"),
         call. = FALSE)
  }
}



# Check for the presence of a given set of labels in an 'ind' file.
# Fail if there a sample was not found.
check_presence <- function(labels, data) {
  not_present <- setdiff(labels, suppressMessages(read_ind(data)$label))
  if (length(not_present) > 0) {
    stop("The following samples are not present in the data: ",
         paste(not_present, collapse = ", "), call. = FALSE)
  }
}



#' Download example SNP data.
#'
#' The data is downloaded to a temporary directory by default.
#' 
#' @param dirname Directory in which to put the data (EIGENSTRAT trio
#'     of snp/geno/ind files).
#'
#' @export
download_data <- function(dirname = tempdir()) {
  dest <- file.path(dirname, "snps.tgz")
  utils::download.file(
    "http://bioinf.eva.mpg.de/admixr/snps.tar.gz",
    destfile = dest,
    method = "wget",
    quiet = TRUE,
    extra = "--no-check-certificate"
  )
  system(paste0("cd ", dirname, "; tar xf ", dest, "; rm snps.tgz"))
  file.path(dirname, "snps", "snps")
}



# this looks incredibly hacky, but seems to be a solution to my R CMD check
# "missing global variable" NOTE woes caused by dplyr code (the following
# are not actually global variables)
utils::globalVariables(
  names = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "chrom", "pos", "snp_id", "ref", "alt", "gen_dist",
            "sample_id", "name", "target", ".", "start", "end",
            "model", "noutgroups", "outgroups", "pattern", "pvalue"),
            package = "admixr")



#' Pipe operator
#' 
#' Added via usethis::use_pipe().
#'
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
