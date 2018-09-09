eigenstrat <- function(prefix) {
    prefix <- path.expand(prefix)
    data <- list(
        ind = paste0(prefix, ".ind"),
        snp = paste0(prefix, ".snp"),
        geno = paste0(prefix, ".geno"),
        group = NULL,
        exclude = NULL
    )
    class(data) <- "eigenstrat"
    data
}

print.eigenstrat <- function(.data) {
    cat(paste(
        "EIGENSTRAT data\n===============",
        "\nind path:", .data$ind, 
        "\nsnp path:", .data$snp, 
        "\ngeno path:", .data$geno,
        "\n\nmodifiers:",
        "\nlabels:", ifelse(is.null(.data$group), "none", .data$group),
        "\nexclude:", ifelse(is.null(.data$exclude), "none", .data$exclude),
        "\n"
    ))
}

process_filter <- function(data, exclude, outfile) {
    # read a table of previously excluded sites, if present
    if (!is.null(data$exclude)) {
        prev_exclude <- read_snp(data$exclude) %>% data.table::setDT()
    } else {
        prev_exclude <- NULL
    }

    # generate a combined set of sites to exclude
    exclude <- rbind(exclude, prev_exclude) %>% unique(., by = c("chrom", "pos"))
    if (nrow(exclude) < nrow(read_snp(data$snp))) {
        write_snp(exclude, outfile)
        data$exclude <- path.expand(outfile)
    } else {
        stop("No sites remaining after the filtering.", call. = FALSE)
    }
    data
}


# read_snp(eigenstrat("all_bigyri_YRI")$snp) %>% sample_n(15) %>% arrange(as.integer(chrom), pos) %>% mutate(start = pos - 1, end = pos) %>% select(chrom, start, end) %>% write_tsv("sites.bed", col_names=F)
filter_bed <- function(data, bed, remove = FALSE, outfile = tempfile()) {
    # process BED file
    bed <- data.table::fread(
        bed,
        col.names = c("chrom", "start", "end"),
        colClasses = c("character", "integer", "integer")
    )
    data.table::setkey(bed, chrom, start, end)

    # process SNP file from the prefix
    snp <- read_snp(data$snp) %>% dplyr::mutate(start = pos - 1, end = pos) %>% data.table::setDT()

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

filter_damage <- function(data, outfile = tempfile()) {
    filtered_snp <- read_snp(data$snp) %>%
        dplyr::filter(
            (ref == "C" & alt == "T") |
            (ref == "T" & alt == "C") |
            (ref == "G" & alt == "A") |
            (ref == "A" & alt == "G")
        )
    data <- process_filter(data, filtered_snp, outfile)
    data
}



#' Change labels of several populations or samples
#'
#' Create a modified EIGENSTRAT ind file with changed population labels.
#'
#' @param data EIGENSTRAT trio.
#' @param labels A named list of labels to merge.
#' @param outfile Where to write the modified ind file.
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



#' Reset modifications to an EIGENSTRAT data
#'
reset <- function(data) {
    data$group <- NULL
    data$exclude <- NULL
    data
}

