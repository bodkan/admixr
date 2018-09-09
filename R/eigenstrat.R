#' EIGENSTRAT data constructor
#'
#' Construct an EIGENSTRAT S3 object from a specified path to the ind/snp/geno
#' trio.
#'
#' The EIGENSTRAT data S3 object encapsulates all paths to data files required
#' for an analysis.
#'
#' @param prefix Shared path to an EIGENSTRAT trio (set of ind/snp/geno files).
#'
#' @return S3 object of a class EIGENSTRAT
#'
#' @export
eigenstrat <- function(prefix) {
    prefix <- path.expand(prefix)
    data <- list(
        ind = paste0(prefix, ".ind"),
        snp = paste0(prefix, ".snp"),
        geno = paste0(prefix, ".geno"),
        group = NULL,
        exclude = NULL
    )
    class(data) <- "EIGENSTRAT"
    data
}



#' EIGENSTRAT print method
#'
#' Print EIGENSTRAT object components.
#'
#' @param data EIGENSTRAT data object.
#'
#' @export
print.EIGENSTRAT <- function(data) {
    cat(paste(
        "EIGENSTRAT data\n===============",
        "\nind path:", data$ind, 
        "\nsnp path:", data$snp, 
        "\ngeno path:", data$geno,
        "\n\nmodifiers:",
        "\nlabels:", ifelse(is.null(data$group), "none", data$group),
        "\nexclude:", ifelse(is.null(data$exclude), "none", data$exclude),
        "\n"
    ))
}



#' Filter EIGENSTRAT data based on a given BED file
#'
#' Keep (or discard) SNPs that overlap (or lie outside of) regions in a given
#' BED file.
#'
#' @param data EIGENSTRAT data object.
#' @param bed Path to a BED file.
#' @param remove Remove sites falling inside the BED file regions? By default,
#'     sites that do not overlap BED regions are removed.
#' @param outfile Path to an output snp file with coordinates of excluded sites.
#'
#' @return Updated S3 EIGENSTRAT data object.
#'
#' @export
filter_bed <- function(data, bed, remove = FALSE, outfile = tempfile()) {
    # process BED file
    bed <- data.table::fread(bed, col.names = c("chrom", "start", "end"),
                             colClasses = c("character", "integer", "integer"))
    data.table::setkey(bed, chrom, start, end)

    # process SNP file from the prefix
    snp <- read_snp(data) %>% dplyr::mutate(start = pos - 1, end = pos) %>% data.table::setDT()

    # get data.table indices of SNPs within/outside given BED regions
    overlap <- data.table::foverlaps(snp, bed, which = TRUE)
    # filter the result based on whether an overlap or a complement is needed
    if (remove) {
        overlap <- overlap[!is.na(overlap$yid), ]
    } else {
        overlap <- overlap[is.na(overlap$yid), ]
    }

    # extract only those sites passing the filter
    site_idx <- unique(overlap$xid)
    exclude <- snp[site_idx, ] %>% dplyr::select(-c(start, end))

    # modify the eigenstrat object and return it
    data <- process_filter(data, exclude, outfile)
    data
}



#' Filter out transitions (C->T and G->A substitutions)
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
#' @export
remove_transitions <- function(data, outfile = tempfile()) {
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
#' @param labels A named list of labels to merge (each new group defined as a
#'     character vector).
#' @param outfile Path to an output snp file with coordinates of excluded sites.
#'
#' @return Updated S3 EIGENSTRAT data object with an additional 'group' slot
#'     specifying the path to a new ind file that will be used in downstream
#'     analysis.
#'
#' @export
relabel <- function(data, labels, outfile = tempfile()) {
  new_lines <- lines <- ifelse(!is.null(data$group), data$group, data$ind) %>% readLines

  # iterate over the lines in the "ind" file, replacing population
  # labels with their substitutes
  for (label in names(labels)) {
    regex <- paste0("(", paste(labels[[label]], collapse = "|"), ")$")
    new_lines <- stringr::str_replace(new_lines, regex, label)
  }

  if (!all(new_lines == lines)) {
      writeLines(new_lines, outfile)
      data$group <- outfile
  } else {
     warning("No labels have been changed")
  }
  data
}



#' Reset modifications to an EIGENSTRAT object
#'
#' Set 'exclude' and 'group' modifications of snp and ind files, respectively,
#' to NULL.
#'
#' @param data EIGENSTRAT data object.
#' @return EIGENSTRAT data S3 object.
#'
#' @export
reset <- function(data) {
    data$group <- NULL
    data$exclude <- NULL
    data
}



# Process the result of SNP filtering and return an updated S3 object
process_filter <- function(data, exclude, outfile) {
    # read a table of previously excluded sites, if present
    if (!is.null(data$exclude)) {
        prev_exclude <- read_snp(data, exclude = TRUE) %>% data.table::setDT()
    } else {
        prev_exclude <- NULL
    }

    # generate a combined set of sites to exclude
    exclude <- rbind(exclude, prev_exclude) %>% unique(., by = c("chrom", "pos"))
    if (nrow(exclude) < nrow(read_snp(data))) {
        write_snp(exclude, outfile)
        data$exclude <- path.expand(outfile)
    } else {
        stop("No sites remaining after the filtering.", call. = FALSE)
    }
    data
}

