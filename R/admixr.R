# Reading output log files --------------------------------------------------


#' Read output log file from a qpF4ratio run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with the parsed results.
#' @export
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
        paste0("\n", collapse="\n") %>%
        readr::read_delim(delim=" ", col_names=FALSE) %>%
        setNames(c("A", "O", "X", "C", "A", "O", "B", "C", "alpha", "stderr", "Zscore")) %>%
        .[c("A", "B", "X", "C", "O", "alpha", "stderr", "Zscore")]

    res_df
}


#' Read output log file from a qpDstat run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with the parsed results.
#' @export
read_qpDstat <- function(file) {
    log_lines <- readLines(file) %>% .[!stringr::str_detect(., "warning")]

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

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        readr::read_delim(delim=" ", col_names=FALSE) %>%
        .[c(1:6, ncol(.) - 2, ncol(.) - 1, ncol(.))] %>% # remove column with "best" if present
        setNames(c("W", "X", "Y", "Z", result_col, "Zscore", "BABA", "ABBA", "n_snps"))

    res_df
}


#' Read output log file from a qp3Pop run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with parsed results.
#' @export
read_qp3Pop <- function(file) {
    log_lines <- readLines(file)

    # extract the number of analyzed population quadruples
    n_quads <- length(log_lines) - (which(stringr::str_detect(log_lines, "^nrows, ncols:"))) - 1

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[stringr::str_detect(log_lines, "result:")] %>%
        stringr::str_replace("result: ", "") %>%
        stringr::str_replace_all(" +", " ") %>%
        stringr::str_replace_all("^ | $", "")

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        readr::read_delim(delim=" ", col_names=FALSE) %>%
        setNames(c("A", "B", "C", "f3", "stderr", "Zscore", "n_snps"))

    res_df
}


# EIGENSTRAT manipulation utilities --------------------------------------------------


#' Merge populations from an EIGENSTRAT "ind" file under a single
#' population label.
#'
#' @param file EIGENSTRAT ind file to modify.
#' @param modified_file Modified EIGENSTRAT ind filename.
#' @param merge List of labels to merge. List names specified labels
#'     to merge into.
#' @export
merge_pops <- function(file, modified_file, merge) {
    # merge=list(ancient_NearEast=merge_what, present_NearEast=c("Yemenite_Jew", "Jordan", "Samaritan", "Bedouin", "Palestinian"))
    lines <- readLines(file)

    # iterate over the lines in the "ind" file, replacing population
    # labels with their substitutes
    for (merge_into in names(merge)) {
        regex <- paste0("(", paste(merge[[merge_into]], collapse="|"), ")$")
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
#' @export
#'
#' @import dplyr
snps_present <- function(geno, prop=FALSE) {
    fn <- ifelse(prop, mean, sum)
    summarise_all(geno, funs(fn(. != 9)))
}


#' Calculate the number (or proportion) of sites with an allele
#' missing for each sample.
#'
#' @param geno EIGENSTRAT geno dataframe.
#' @param prop Calculate the proportion of missing alleles
#'     instead.
#'
#' @return A named vector of counts or proportions.
#' @export
#'
#' @import dplyr
snps_missing <- function(geno, prop=FALSE) {
    fn <- ifelse(prop, mean, sum)
    summarise_all(geno, funs(fn(. == 9)))
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
#' @import dplyr
#' @export
subset_sites <- function(prefix, out_prefix, bed_file, complement=FALSE) {
    coords <- readr::read_table2(bed_file, col_names=c("chrom", "start", "pos"),
                          col_types="cii", progress=FALSE) %>% select(-start)

    geno <- read_geno(paste0(prefix, ".geno"))
    snp <- read_snp(paste0(prefix, ".snp"))
    combined <- bind_cols(snp, geno)

    # determine which function to call on the coordinates
    fun <- ifelse(complement, anti_join, inner_join)
    combined_subset <- fun(combined, coords, by=c("chrom", "pos"))

    # write the new snp file
    write_tsv(select(combined_subset, id:alt), path=paste0(out_prefix, ".snp"), col_names=FALSE)
    # write the new geno file
    writeLines(apply(select(combined_subset, -(id:alt)), 1, paste, collapse=""),
               con=paste0(out_prefix, ".geno"))
    # write the new ind file
    invisible(file.copy(from=paste0(prefix, ".ind"), to=paste0(out_prefix, ".ind")))
}


#' Convert EIGENSTRAT files into a VCF file.
#'
#' This function reads the genotypes from a 'geno' file, their
#' coordinates and reference/alternative alleles and the individuals'
#' names from an 'ind' file and generates a minimalistic VCF file.
#'
#' Compressing and indexing requires having "bgzip" and "tabix"
#' commands in $PATH.
#'
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param vcf_file Path to the VCF file that will be generated.
#' @param compress Compress the VCF with bgzip?
#' @param index Index the VCF with tabix?
#'
#' @export
eigenstrat_to_vcf <- function(prefix, vcf_file, compress=TRUE, index=TRUE) {
    geno <- read_geno(paste0(prefix, ".geno"), paste0(prefix, ".ind"))
    snp <- read_snp(paste0(prefix, ".snp"))
    ind <- read_ind(paste0(prefix, ".ind"))

    # construct a minimal VCF header
    header <- c("##fileformat=VCFv4.1",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                sapply(unique(snp$chrom), function(chrom) {
                    paste0("##contig=<ID=", chrom, ">")
                }))

    # generate a dataframe with the "body" of the VCF file (info and GT columns)
    info_cols <- mutate(snp, ID=".", QUAL="0", FILTER=".", INFO=".", FORMAT="GT") %>%
        select(`#CHROM`=chrom, POS=pos, ID, REF=ref, ALT=alt, QUAL, FILTER, INFO, FORMAT)
    gt_cols <- mutate_all(geno, eigenstrat_to_gt)
    body_cols <- bind_cols(info_cols, gt_cols)

    writeLines(header, vcf_file)
    write_tsv(body_cols, vcf_file, col_names=TRUE, append=TRUE)

    if (compress) system(paste("bgzip", vcf_file))
    if (index)    system(paste("tabix", paste0(vcf_file, ".gz")))
}


#' Convert VCF file into three-file EIGENSTRAT format.
#'
#' @param vcf_file Path to the VCF file.
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path) that will be generated.
#'
#' @export
vcf_to_eigenstrat <- function(vcf_file, prefix) {
    vcf <- readr::read_tsv(vcf_file, comment="##") %>%
        rename(chrom=`#CHROM`, pos=POS, ref=REF, alt=ALT) %>%
        select(-c(ID, QUAL, FILTER, INFO, FORMAT))

    # generate dataframes with the 3 EIGENSTRAT info tables
    snp <- select(vcf, chrom, pos, ref, alt) %>%
        mutate(snp_id=paste(chrom, pos, sep="_"), gen_dist="0.0") %>%
        select(snp_id, chrom, gen_dist, pos, ref, alt)
    ind <- tibble(
        sample_id=select(vcf, -c(chrom, pos, ref, alt)) %>% names,
        sex="U",
        label=sample_id
    )
    geno <- select(vcf, -c(chrom, pos, ref, alt)) %>%
        mutate_all(gt_to_eigenstrat)

    # write all three EIGENSTRAT files
    write_tsv(snp, paste0(prefix, ".snp"), col_names=FALSE)
    write_tsv(ind, paste0(prefix, ".ind"), col_names=FALSE)
    writeLines(apply(geno, 1, paste, collapse=""), paste0(prefix, ".geno"))
}
