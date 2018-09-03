#' Merge multiple samples/populations under a single label
#'
#' @param ind Path to an original ind file.
#' @param modified_ind Path to a final modified ind file.
#' @param labels A named list of labels to merge.
#'
#' @examples
#'
#' \dontrun{
#' # This will create a new ind file with labels in the 3rd column replaced by
#' # Europe", "EastAfrica" and "WestAfrica", respectively.
#' group_labels(
#'     ind = IND_FILE,
#'     modified_ind = paste0(IND_FILE, ".merged"),
#'     labels = list(Europe = c("French", "Sardinian", "Czech"),
#'                   WestAfrica = c("Yoruba", "Mende"))
#' )
#' }
#'
#' @export
group_labels <- function(ind, modified_ind, labels) {
  lines <- readLines(ind)

  # iterate over the lines in the "ind" file, replacing population
  # labels with their substitutes
  for (merge_into in names(labels)) {
    regex <- paste0("(", paste(labels[[merge_into]], collapse = "|"), ")$")
    lines <- stringr::str_replace(lines, regex, merge_into)
  }

  writeLines(lines, modified_ind)
}


#' Merge two sets of EIGENSTRAT datasets
#' 
#' This function utilizes the 'mergeit' command from ADMIXTOOLS.
#'
#' @param prefix A prefix of a final merged EIGENSTRAT dataset.
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


#' Count the number/proportion of present/missing sites in each sample
#'
#' @param prefix An EIGENSTRAT prefix.
#' @param prop Calculate the proportion instead of counts?
#' @param missing Count present SNPs or missing SNPs?
#'
#' @return A data.frame object with SNP counts/proportions.
#'
#' @export
#' @import rlang
count_snps <- function(prefix, missing = FALSE, prop = FALSE) {
    fn <- ifelse(prop, mean, sum)
    if (missing) {
        op <- `==`
        col <- "missing"
    } else {
        op <- `!=`
        col <- "present"
    }
    eigenstrat <- read_eigenstrat(prefix)
    dplyr::summarise_all(eigenstrat$geno, dplyr::funs(fn(op(., 9)))) %>%
        tidyr::gather(name, !!col)
}



#' Generate coordinates of SNPs that overlap/exclude regions in BED file
#'
#' @param snp Path to a snp file.
#' @param bed Path to a BED file.
#' @param output Path to an output snp file with coordinates to exclude.
#' @param exclude Exclude sites falling inside the BED file regions?
#'
#' @export
filter_sites <- function(snp, bed, output, exclude = FALSE) {
  # process BED file
  dt_bed <- data.table::fread(
    bed,
    col.names = c("chrom", "start", "end"),
    colClasses = c("character", "integer", "integer")
  )
  data.table::setkey(dt_bed, chrom, start, end)

  # process SNP file from the prefix
  dt_snp <- read_snp(snp) %>%
    dplyr::mutate(chrom = as.character(chrom), start = pos - 1, end = pos) %>%
    data.table::setDT()
  data.table::setkey(dt_snp, chrom, start, end)                                          
                                                                                      
  # get data.table indices of SNPs within/outside given BED regions                   
  overlap <- data.table::foverlaps(dt_snp, dt_bed, which = TRUE)                            
  # filter the result based on whether an overlap or a complement is needed           
  if (exclude)                                                                        
    overlap <- overlap[is.na(overlap$yid), ]                                         
  else
    overlap <- overlap[!is.na(overlap$yid), ]

  # extract only those sites passing the filter
  site_idx <- unique(overlap$xid)
  if (!length(site_idx)) stop("No sites remaining after the overlap!")
  snp_subset <- dt_snp[site_idx, ]

  dplyr::select(snp_subset, -c(start, end)) %>% write_snp(output)
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

#' Download example SNP data.
#'
#' @param dirname Directory in which to put the data (EIGENSTRAT trio of snp/geno/ind files).
#'
#' @export
download_data <- function(dirname = tempdir()) {
    dest <- file.path(dirname, "snps.tgz")
    utils::download.file(
        "https://www.dropbox.com/s/bw1eswcp2domisv/snps.tgz?dl=0",
        destfile = dest,
        method = "wget",
        quiet = TRUE
    )
    system(paste0("cd ", dirname, "; tar xf ", dest, "; rm snps.tgz"))
    file.path(dirname, "snps", "snps")
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


# this looks incredibly hacky, but seems to be a solution to my R CMD check
# "missing global variable" NOTE woes caused by dplyr code (the following
# are not actually global variables)
utils::globalVariables(
    names = c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
              "FORMAT", "chrom", "pos", "snp_id", "ref", "alt", "gen_dist",
              "sample_id", "name", "target", ".", "start", "end"),
    package = "admixr"
)
