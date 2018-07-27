#' Merge populations from an EIGENSTRAT "ind" file under a single
#' population label.
#'
#' @param ind EIGENSTRAT ind file to modify.
#' @param modified_ind Modified EIGENSTRAT ind filename.
#' @param pop_labels Named list of labels to merge. Names specify labels
#'     to merge into.
#'
#' @examples
#'
#' # This will create a new ind file with labels in the 3rd column replaced by
#' # "Europe", "EastAfrica" and "WestAfrica", respectively.
#' \dontrun{merge_labels(
#'     ind = IND_FILE,
#'     modified_ind = paste0(IND_FILE, ".merged"),
#'     pop_labels = list(Europe = c("French", "Sardinian", "Czech",
#'                       EastAfrica = "Dinka",
#'                       WestAfrica = c("Mbuti", "Mende))
#' )}
#'
#' @export
merge_labels <- function(ind, modified_ind, pop_labels) {
  lines <- readLines(ind)

  # iterate over the lines in the "ind" file, replacing population
  # labels with their substitutes
  for (merge_into in names(pop_labels)) {
    regex <- paste0("(", paste(merge[[merge_into]], collapse = "|"), ")$")
    lines <- stringr::str_replace(lines, regex, merge_into)
  }

  writeLines(lines, modified_file)
}



# Filtering functions  --------------------------------------------------


#' Calculate the number (or proportion) of sites with an allele
#' present (i.e. not 9) for each sample.
#'
#' @param geno EIGENSTRAT geno dataframe.
#' @param prop Calculate the proportion of non-missing alleles
#'     instead.
#'
#' @return A named vector of counts or proportions.
#'
#' @export
snps_present <- function(geno, prop = FALSE) {
  fn <- ifelse(prop, mean, sum)
  dplyr::summarise_all(geno, funs(fn(. != 9)))
}


#' Calculate the number (or proportion) of sites with an allele
#' missing for each sample.
#'
#' @param geno EIGENSTRAT geno dataframe.
#' @param prop Calculate the proportion of missing alleles
#'     instead.
#'
#' @return A named vector of counts or proportions.
#'
#' @export
snps_missing <- function(geno, prop = FALSE) {
  fn <- ifelse(prop, mean, sum)
  dplyr::summarise_all(geno, funs(fn(. == 9)))
}


#' Create a new set of EIGENSTRAT files by intersecting the original
#' data with a given set of coordinates.
#'
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param out_prefix Prefix of the generated EIGENSTRAT files with the
#'     subset of the data.
#' @param bed_file Path to the 3 column BED file to intersect with.
#' @param complement Perform an intersect or a complement operation?
#'
#' @export
#'
#' @importFrom magrittr "%>%"
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
    if (!cmd %in% c("qpDstat", "qpF4ratio", "qp3Pop", "qpAdm")) {
        stop("ADMIXTOOLS command '", cmd, "' is not supported or does not exist")
    }

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
