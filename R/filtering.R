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
    # process BED file
    dt_bed <- data.table::fread(
        bed,
        col.names = c("chrom", "start", "end"),
        colClasses = c("character", "integer", "integer")
    )
    data.table::setkey(dt_bed, chrom, start, end)

    # process SNP file from the prefix
    snp <- paste0(prefix, ".snp")
    dt_snp <- read_snp(snp) %>%
        dplyr::mutate(chrom = as.character(chrom), start = pos - 1, end = pos) %>%
        data.table::setDT()
    data.table::setkey(dt_snp, chrom, start, end)                                          

    # get data.table indices of SNPs within/outside given BED regions                   
    overlap <- data.table::foverlaps(dt_snp, dt_bed, which = TRUE)                            
    # filter the result based on whether an overlap or a complement is needed           
    if (remove)                                                                        
        overlap <- overlap[is.na(overlap$yid), ]                                         
    else
        overlap <- overlap[!is.na(overlap$yid), ]

    # extract only those sites passing the filter
    site_idx <- unique(overlap$xid)
    if (!length(site_idx)) stop("No sites remaining after the overlap!")
    snp_subset <- dt_snp[site_idx, ]

    dplyr::select(snp_subset, -c(start, end)) %>% write_snp(outsnp)
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
filter_damage <- function(prefix, outsnp) {
    read_snp(paste0(prefix, ".snp")) %>%
        dplyr::filter(
            (ref == "C" & alt == "T") |
            (ref == "T" & alt == "C") |
            (ref == "G" & alt == "A") |
            (ref == "A" & alt == "G")
        ) %>%
        write_snp(outsnp)
}

