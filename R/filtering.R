#' Generate coordinates of SNPs that overlap/exclude regions in BED file
#'
#' This function normally produces a snp file with sites that did not overlap
#' regions in the specified BED file, that can be then directly supplied to
#' any admixr statistic function using it's 'exclude = ' argument. If the
#' BED file contains regions that have to be filtered out, set 'remove = TRUE'.
#'
#' @param prefix An EIGENSTRAT prefix.
#' @param bed Path to a BED file.
#' @param outsnp Path to an output snp file with coordinates of excluded sites.
#' @param remove Remove sites falling inside the BED file regions? By default, sites
#'     that do not overlap BED regions are removed.
#'
#' @export
filter_sites <- function(prefix, bed, outsnp, remove = FALSE) {


}



#' Generate a snp file with positions of transversion SNPs.
#'
#' The output of this function can be used as 'exclude = ' argument in all main
#' admixr statistics functions.
#'
#' @param prefix An EIGENSTRAT prefix.
#' @param outsnp Path to an output snp file with coordinates of excluded sites.
#'
#' @export

