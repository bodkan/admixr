#' Merge multiple samples/populations under a single label.
#'
#' @param ind EIGENSTRAT ind file to modify.
#' @param modified_ind Modified EIGENSTRAT ind filename.
#' @param labels Named list of labels to merge. Names specify labels
#'     to merge into.
#'
#' @examples
#'
#' \dontrun{
#' # This will create a new ind file with labels in the 3rd column replaced by
#' # Europe", "EastAfrica" and "WestAfrica", respectively.
#' merge_labels(
#'     ind = IND_FILE,
#'     modified_ind = paste0(IND_FILE, ".merged"),
#'     labels = list(Europe = c("French", "Sardinian", "Czech"),
#'                       WestAfrica = c("Yoruba", "Mende"))
#' )
#' }
#'
#' @export
merge_labels <- function(ind, modified_ind, labels) {
  lines <- readLines(ind)

  # iterate over the lines in the "ind" file, replacing population
  # labels with their substitutes
  for (merge_into in names(labels)) {
    regex <- paste0("(", paste(labels[[merge_into]], collapse = "|"), ")$")
    lines <- stringr::str_replace(lines, regex, merge_into)
  }

  writeLines(lines, modified_ind)
}


#' Merge two sets of EIGENSTRAT datasets (utilizing the 'mergeit' command).
#'
#' @param prefix Prefix of the merged dataset.
#' @param input1,input2 Prefixes of two EIGENSTRAT datasets to merge.
#'
#' @export
merge_eigenstrat <- function(prefix, input1, input2) {
  parfile <- tempfile()
  paste0(
    "outputformat: EIGENSTRAT\n",
    "strandcheck: NO\n",
    "geno1: ", input1, ".geno\n",
    "snp1: ", input1, ".snp\n",
    "ind1: ", input1, ".ind\n",
    "geno2: ", input2, ".geno\n",
    "snp2: ", input2, ".snp\n",
    "ind2: ", input2, ".ind\n",
    "genooutfilename: ", prefix, ".geno\n",
    "snpoutfilename: ", prefix, ".snp\n",
    "indoutfilename: ", prefix, ".ind"
  ) %>% writeLines(text = ., con = parfile)

  return_value <- run_cmd("mergeit", parfile, "/dev/null")
  if (return_value) cat("\nMerge command ended with an error -- see above.\n")
}

# Filtering functions  --------------------------------------------------


#' Calculate the number/proportion of covered sites in each sample.
#'
#' @param prefix EIGENSTRAT prefix.
#' @param prop Calculate the proportion of non-missing alleles instead?
#'
#' @return data.frame object with SNP counts/proportions.
#'
#' @export
snps_present <- function(prefix, prop = FALSE) {
  fn <- ifelse(prop, mean, sum)
  eigenstrat <- read_eigenstrat(prefix)
  dplyr::summarise_all(eigenstrat$geno, dplyr::funs(fn(. != 9))) %>%
    tidyr::gather(name, nsnps)
}


#' Calculate the number/proportion of missing sites in each sample.
#'
#' @param prefix EIGENSTRAT prefix.
#' @param prop Calculate the proportion of non-missing alleles instead?
#'
#' @return data.frame object with SNP counts/proportions.
#'
#' @export
snps_missing <- function(prefix, prop = FALSE) {
  fn <- ifelse(prop, mean, sum)
  eigenstrat <- read_eigenstrat(prefix)
  dplyr::summarise_all(eigenstrat$geno, dplyr::funs(fn(. == 9))) %>%
    tidyr::gather(name, nsnps)
}


#' Generate new EIGENSTRAT dataset overlapping a given BED file.
#'
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param out_prefix Prefix of the generated EIGENSTRAT files with the
#'     subset of the data.
#' @param bed_file Path to the 3 column BED file to intersect with.
#' @param complement Perform an intersect or a complement operation?
#'
#' @export
subset_sites <- function(prefix, out_prefix, bed_file, complement = FALSE) {
  coords <- readr::read_table2(
    bed_file,
    col_names = c("chrom", "start", "pos"),
    col_types = "cii",
    progress = FALSE
  ) %>%
    dplyr::select(-start)
  
  geno <- read_geno(paste0(prefix, ".geno"))
  snp <- read_snp(paste0(prefix, ".snp"))
  combined <- dplyr::bind_cols(snp, geno)
  
  # determine which function to call on the coordinates
  fun <- ifelse(complement, dplyr::anti_join, dplyr::inner_join)
  combined_subset <- fun(combined, coords, by = c("chrom", "pos"))
  
  # write the new snp file
  dplyr::select(combined_subset, id:alt) %>%  
    readr::write_tsv(path = paste0(out_prefix, ".snp"), col_names = FALSE)
  # write the new geno file
  dplyr::select(combined_subset, -(id:alt)) %>%
    apply(1, paste, collapse = "") %>%
    writeLines(con = paste0(out_prefix, ".geno"))
  # write the new ind file
  invisible(file.copy(from = paste0(prefix, ".ind"),
                      to = paste0(out_prefix, ".ind")))
}


# Run a specified ADMIXTOOLS command.
run_cmd <- function(cmd, par_file, log_file) {
  system(paste(cmd, "-p", par_file, ">", log_file))
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


# Check for the presence of a given set of labels in an 'ind' file.
# Fail if there a sample was not found.
check_presence <- function(labels, prefix = NULL, ind = NULL) {
    if (!is.null(ind)) {
        path <- ind
    } else {
        path <- paste0(prefix, ".ind")
    }

    not_present <- setdiff(labels, suppressMessages(read_ind(path)$label))
    if (length(not_present) > 0) {
        stop("The following samples are not present in '", ind, "': ",
             paste(not_present, collapse = ", "))
    }
}


# Look for path to the ADMIXTOOLS directory.
admixtools_path <- function() {
  system("which qpDstat", intern = TRUE) %>% stringr::str_replace("/bin.*", "")
}

#' Pipe operator
#' 
#' Added via usethis::use_pipe().
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
