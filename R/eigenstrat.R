#' EIGENSTRAT data constructor
#'
#' This function creates an instance of the EIGENSTRAT S3 class, which
#' encapsulates all paths to data files required for an ADMIXTOOLS analysis.
#'
#' @param prefix Shared path to an EIGENSTRAT trio (set of ind/snp/geno files).
#' @param ind,snp,geno Paths to individual EIGENSTRAT components.
#' @param exclude Pre-defined snp file with excluded sites.
#'
#' @return S3 object of the EIGENSTRAT class.
#'
#' @examples
#' \dontrun{# download an example genomic data and get the path prefix to the
#' # trio of snp/geno/ind files in an EIGENSTRAT format
#' prefix <- download_data(dirname = tempdir())
#'
#' # wrap the trio of snp/geno/ind files in an object of the class
#' # EIGENSTRAT
#' snps <- eigenstrat(prefix)
#' }
#'
#' @export
eigenstrat <- function(prefix = NULL, ind = NULL, snp = NULL, geno = NULL, exclude = NULL) {
  if (is.null(prefix) & any(is.null(c(ind, snp, geno))))
    stop("Insufficient information to get paths to ind/snp/geno files")

  data <- list(
    ind = path.expand(ifelse(!is.null(ind), ind, paste0(prefix, ".ind"))),
    snp = path.expand(ifelse(!is.null(snp), snp, paste0(prefix, ".snp"))),
    geno = path.expand(ifelse(!is.null(geno), geno, paste0(prefix, ".geno"))),
    group = NULL,
    exclude = NULL
  )

  if (!all(sapply(c(data$ind, data$snp, data$geno), file.exists)))
    stop("Not all three ind/snp/geno files present", call. = FALSE)

  if (!is.null(exclude)) {
    data <- add_excluded_snps(data, exclude)
  }

  class(data) <- "EIGENSTRAT"
  data
}



#' EIGENSTRAT print method
#'
#' Print EIGENSTRAT object components.
#'
#' @param x EIGENSTRAT data object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
print.EIGENSTRAT <- function(x, ...) {
  cat(paste0(
    "EIGENSTRAT object\n",
    "=================\n",
    "components:",
    "\n  ind file: ", x$ind,
    "\n  snp file: ", x$snp,
    "\n  geno file: ", x$geno,
    "\n"))

  if (any(!is.null(c(x$group, x$exclude)))) cat("\nmodifiers:\n")

  if (!is.null(x$group)) cat("  groups: ", x[["group"]], "\n")

  if (!is.null(x$exclude)) {
    cat("  excluded sites: ", x$exclude, "\n")
    cat("      (SNPs excluded: ", x$n_excluded,
             ", SNPs remaining: ", x$n_included, ")\n", sep = "")
  }
}



#' Merge two sets of EIGENSTRAT datasets
#' 
#' This function utilizes the 'mergeit' command distributed in ADMIXTOOLS.
#'
#' @param merged Prefix of the path to the merged EIGENSTRAT snp/ind/geno trio.
#' @param a,b Two EIGENSTRAT objects to merge.
#' @param strandcheck Deal with potential strand issues? Mostly for historic reasons. For details see the README of ADMIXTOOLS convertf.
#'
#' @export
merge_eigenstrat <- function(merged, a, b, strandcheck = "NO") {
  check_type(a, "EIGENSTRAT")
  check_type(b, "EIGENSTRAT")

  parfile <- tempfile()
  paste0(
    "outputformat: EIGENSTRAT\n",
    "strandcheck: ", strandcheck, "\n",
    "geno1: ", a$geno, "\n",
    "snp1: ", a$snp, "\n",
    "ind1: ", a$ind, "\n",
    "geno2: ", b$geno, "\n",
    "snp2: ", b$snp, "\n",
    "ind2: ", b$ind, "\n",
    "genooutfilename: ", merged, ".geno\n",
    "snpoutfilename: ", merged, ".snp\n",
    "indoutfilename: ", merged, ".ind"
  ) %>% writeLines(text = ., con = parfile)

  return_value <- run_cmd("mergeit", parfile, "/dev/null")
  if (return_value) stop("\nMerge command ended with an error -- see above.\n")

  eigenstrat(merged)
}



#' Filter EIGENSTRAT data based on a given BED file
#'
#' Keep (or discard) SNPs that overlap (or lie outside of) regions in a given
#' BED file.
#'
#' This function requires a functioning bedtools installation! See:
#'
#' - https://github.com/arq5x/bedtools2
#'
#' - https://bedtools.readthedocs.io/
#'
#' @param data EIGENSTRAT data object.
#' @param bed Path to a BED file.
#' @param remove Remove sites falling inside the BED file regions? By default,
#'     sites that do not overlap BED regions are removed.
#' @param outfile Path to an output snp file with coordinates of excluded sites.
#' @param bedtools_args Optional arguments to `bedtools intersect` such as \code{"-sorted"}
#'   or \code{"-sorted -nonamecheck"}.
#'
#' @return Updated S3 EIGENSTRAT data object.
#'
#' @examples
#' \dontrun{# download an example genomic data set
#' prefix <- download_data(dirname = tempdir())
#' # create an EIGENSTRAT R object from the downloaded data
#' snps <- eigenstrat(prefix)
#'
#' # get the path to an example BED file
#' bed <- file.path(dirname(prefix), "regions.bed")
#'
#' # BED file contains regions to keep in an analysis
#' snps_kept <- filter_bed(snps, bed)
#' # BED file contains regions to remove from an analysis
#' snps_removed <- filter_bed(snps, bed, remove = TRUE)
#' }
#'
#' @export
filter_bed <- function(data, bed, remove = FALSE, outfile = tempfile(fileext = ".snp"), bedtools_args = "") {
  check_type(data, "EIGENSTRAT")

  if (file.exists(outfile)) {
    warning(paste("snp file", outfile, "already exists - adding it to the EIGENSTRAT object directly."), call. = FALSE)
    data <- add_excluded_snps(data, outfile)
    return(data)
  }

  if (Sys.which("bedtools") == "")
    stop("bedtools is required for filtering, but is not in your $PATH")

  if (!file.exists(bed))
    stop("BED file '", bed, "' not found", call. = FALSE)

  snp <- read_snp(data)

  tmpbed <- tempfile()
  snp %>%
    dplyr::mutate(start = pos - 1, end = pos) %>%
    dplyr::select(chrom, start, pos) %>%
    readr::write_tsv(tmpbed, col_names = FALSE)

  # get positions of SNPs within/outside given BED regions (the -c
  # argument counts how many times does a site from `a` overlap anything in `b`;
  # for our purposes, this means either 0 (no overlap) or 1 (overlap)
  output <- tempfile()
  # run bedtools
  sprintf("bedtools intersect %s -c -a %s -b %s %s > %s",
          ifelse(remove, "-v", ""), tmpbed, bed, bedtools_args, output) %>% system()
  # collect the results
  snp_hits <- readr::read_tsv(output,
                              col_names = c("chrom", "start", "end", "hit"),
                              col_types = "ciii",
                              progress = FALSE)

  # extract only those sites passing the filter
  exclude <- snp[snp_hits$hit == as.integer(remove), ]

  # modify the eigenstrat object and return it
  data <- process_filter(data, exclude, outfile)
  data
}



#' Remove transversions (C->T and G->A substitutions)
#'
#' Remove substitutions that are more likely to be a result of ancient DNA
#' damage (C->T and G->A substitutions).
#'
#' @param data EIGENSTRAT data object.
#' @param outfile Path to an output snp file with coordinates of excluded sites.
#'
#' @return Updated S3 EIGENSTRAT data object with an additional 'exclude' slot
#'     specifying the path to the set of SNPs to be removed from a downstream
#'     analysis.
#'
#' @examples
#' \dontrun{# download an example genomic data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # perform the calculation only on transversions
#' snps_tv <- transversions_only(snps)
#' results_d <- d(W = "French", X = "Dinka", Y = "Altai", Z = "Chimp", data = snps_tv)
#' }
#'
#' @export
transversions_only <- function(data, outfile = tempfile(fileext = ".snp")) {
  check_type(data, "EIGENSTRAT")

  if (file.exists(outfile)) {
    warning(paste0("snp file '", outfile, "' already exists - adding it to the EIGENSTRAT object directly..."), call. = FALSE)
    data <- add_excluded_snps(data, outfile)
    return(data)
  }

  exclude <- read_snp(data) %>%
    dplyr::filter(
      (ref == "C" & alt == "T") |
        (ref == "T" & alt == "C") |
        (ref == "G" & alt == "A") |
        (ref == "A" & alt == "G")
    )
  data <- process_filter(data, exclude, outfile)
  data
}



#' Change labels of populations or samples
#'
#' Replace population/sample names with specified group labels.
#'
#' @param data EIGENSTRAT trio.
#' @param ... Population/sample names to merge (each new group defined as a
#'     character vector).
#' @param outfile Path to an output snp file with coordinates of excluded sites.
#'
#' @return Updated S3 EIGENSTRAT data object with an additional 'group' slot
#'     specifying the path to a new ind file.
#'#'
#' @examples
#' \dontrun{# download an example genomic data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # group individual samples into larger populations, creating a new
#' # EIGENSTRAT R object
#' new_snps <- relabel(
#'     snps,
#'     European = c("French", "Sardinian"),
#'     African = c("Dinka", "Yoruba", "Mbuti", "Khomani_San"),
#'     Archaic = c("Vindija", "Altai", "Denisova")
#' )
#' }
#'
#' @export
relabel <- function(data, ..., outfile = tempfile(fileext = ".ind")) {
  check_type(data, "EIGENSTRAT")

  labels <- list(...)

  # iterate over the lines in the "ind" file, replacing population
  # labels with their substitutes
  new_lines <- lines <- ifelse(!is.null(data$group), data$group, data$ind) %>% readLines
  for (label in names(labels)) {
    regex <- paste0("(", paste(labels[[label]], collapse = "|"), ")$")
    new_lines <- stringr::str_replace(new_lines, regex, label)
  }

  if (!all(new_lines == lines)) {
    writeLines(new_lines, outfile)
    data$group <- outfile
  } else {
    warning("No labels have been changed", call. = FALSE)
  }
  data
}



#' Reset modifications to an EIGENSTRAT object
#'
#' Set 'exclude' and 'group' modifications of snp and ind files, effectively
#' resetting the dataset into its original state.
#'
#' @param data EIGENSTRAT data object.
#' @return EIGENSTRAT data S3 object.
#'
#'
#' @examples
#' \dontrun{# download an example genomic data set and prepare it for analysis
#' snps <- eigenstrat(download_data(dirname = tempdir()))
#'
#' # group individual samples into larger populations, creating a new
#' # EIGENSTRAT R object
#' new_snps <- relabel(
#'     snps,
#'     European = c("French", "Sardinian"),
#'     African = c("Dinka", "Yoruba", "Mbuti", "Khomani_San"),
#'     Archaic = c("Vindija", "Altai", "Denisova")
#' )
#'
#' # remove the population grouping in the previous step - this
#' # results in the same EIGENSTRAT object tht we started with
#' original_snps <- reset(new_snps)
#' }
#'
#' @export
reset <- function(data) {
  check_type(data, "EIGENSTRAT")

  data$group <- NULL
  data$exclude <- NULL
  data$excluded_n <- NULL
  data$included_n <- NULL
  data
}



# Process the result of SNP filtering and return an updated S3 object
process_filter <- function(data, exclude, outfile) {
  # read a table of previously excluded sites, if present
  if (!is.null(data$exclude)) {
    prev_exclude <- read_snp(data, exclude = TRUE)
  } else {
    prev_exclude <- NULL
  }

  # generate a combined set of sites to exclude
  exclude <- rbind(exclude, prev_exclude) %>% unique(., by = c("chrom", "pos"))
  if (nrow(exclude) < nrow(read_snp(data))) {
    write_snp(exclude, outfile)
    data <- add_excluded_snps(data, outfile)
  } else {
    stop("No sites remaining after the filtering.", call. = FALSE)
  }
  data
}



# Add an 'exclude' field to the EIGENSTRAT object and count remaining
# and excluded SNPs.
add_excluded_snps <- function(data, snpfile) {
    data$exclude <- path.expand(snpfile)
    ## data$n_excluded <- count_lines(data$exclude)
    ## data$n_included <- count_lines(data$snp) - data$n_excluded
    data$n_excluded <- nrow(read_snp(data, exclude = TRUE))
    data$n_included <- nrow(read_snp(data)) - data$n_excluded
    data
}



# Quick hack to count lines in a given file.
# Used just for the purpose of counting the number of SNPs before/after
# filtering.
count_lines <- function(file) {
    wc_out <- system(paste0("wc -l ", file), intern = TRUE)
    as.integer(gsub("^ +([0-9]+) .*", "\\1", wc_out))
}
