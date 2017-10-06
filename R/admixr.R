# High level wrapper functions --------------------------------------------------


#' Run an f4-ratio analysis and return the results as a data.frame.
#'
#' @param X, A, B, C, O Population names, using the terminology of
#'     Patterson et al., 2012
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix' argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @export
qpF4ratio <- function(X, A, B, C, O,
                     prefix=NULL, geno=NULL, snp=NULL, ind=NULL, badsnp=NULL,
                     dir_name=NULL) {
    check_presence(c(X, A, B, C, O), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qpF4ratio__", A, "_", B, "_", C, "_", O)
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qpF4ratio_pop_file(X=X, A=A, B=B, C=C, O=O, file=files[["pop_file"]])
    create_par_file(files[["par_file"]], files[["pop_file"]],
                    prefix, geno, snp, ind, badsnp)

    run_cmd("qpF4ratio", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qpF4ratio(files[["log_file"]]) %>% mutate(setup=setup)
}


#' Calculate D statistics or F4 statistics (which is just the
#' numerator of a D statistic) and return the results as a data.frame.
#'
#' @param W, X, Y, Z Population names, using the terminology of
#'     Patterson et al., 2012
#' @param prefix Prefix of the geno/snp/ind files (including the whole
#'     path).
#' @param geno Path to the genotype file. Overrides the 'prefix' argument.
#' @param snp Path to the snp file. Overrides the 'prefix' argument.
#' @param ind Path to the ind file. Overrides the 'prefix' argument.
#' @param badsnp SNP file with information about ignored sites.
#' @param dir_name Where to put all generated files (temporary
#'     directory by default).
#' @export
qpDstat <- function(W, X, Y, Z,
                    prefix=NULL, geno=NULL, snp=NULL, ind=NULL, badsnp=NULL,
                    dir_name=NULL, f4mode=FALSE) {
    check_presence(c(W, X, Y, Z), prefix, ind)

    # get the path to the population, parameter and log files
    setup <- paste0("qpDstat")
    config_prefix <- paste0(setup, "__", as.integer(runif(1, 0, .Machine$integer.max)))
    files <- get_files(dir_name, config_prefix)

    create_qpDstat_pop_file(W, X, Y, Z, file=files[["pop_file"]])
    create_par_file(files[["par_file"]], files[["pop_file"]],
                    prefix, geno, snp, ind, badsnp, f4mode)

    run_cmd("qpDstat", par_file=files[["par_file"]], log_file=files[["log_file"]])

    read_qpDstat(files[["log_file"]])
}


# Reading output log files --------------------------------------------------


#' Read output log file from a qpF4ratio run.
#'
#' @param file Name of the output log file.
#'
#' @return Tibble object with the parsed results.
#' @export
#'
#' @import stringr readr
read_qpF4ratio <- function(file) {
    log_lines <- readLines(file) %>% .[!str_detect(., "warning")]

    # extract the number of analyzed test populations/individuals
    # (corresponding to the number of rows of the results table)
    n_pops <- log_lines[which(str_detect(log_lines, "^nplist:"))] %>%
        str_extract("[0-9]+$") %>%
        as.integer

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[(length(log_lines) - n_pops) : (length(log_lines) - 1)] %>%
        str_replace("result: ", "") %>%
        str_replace_all(":", "") %>%
        str_replace_all(" +", " ") %>%
        str_replace("^ ", "")

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        read_delim(delim=" ", col_names=FALSE) %>%
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
#'
#' @import stringr readr
read_qpDstat <- function(file) {
    log_lines <- readLines(file) %>% .[!str_detect(., "warning")]

    # extract the number of analyzed population quadruples
    n_quads <- length(log_lines) - (which(str_detect(log_lines, "^nrows, ncols:"))) - 1

    # parse the lines of the results section and extract the names of
    # tested populations/individuals, estimated admixture proportions
    # alpha, std. errors and Z-score
    res_lines <- log_lines[(length(log_lines) - n_quads) : (length(log_lines) - 1)] %>%
        str_replace("result: ", "") %>%
        str_replace_all(" +", " ") %>%
        str_replace_all("^ | $", "")

    result_col <- ifelse(any(str_detect(log_lines, "f4mode: YES")), "f4", "D")

    res_df <- res_lines %>%
        paste0("\n", collapse="\n") %>%
        read_delim(delim=" ", col_names=FALSE) %>%
        .[c(1:6, ncol(.) - 2, ncol(.) - 1, ncol(.))] %>% # remove column with "best" if present
        setNames(c("W", "X", "Y", "Z", result_col, "Zscore", "BABA", "ABBA", "n_snps"))

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
#'
#' @import stringr
merge_pops <- function(file, modified_file, merge) {
    # merge=list(ancient_NearEast=merge_what, present_NearEast=c("Yemenite_Jew", "Jordan", "Samaritan", "Bedouin", "Palestinian"))
    lines <- readLines(file)

    # iterate over the lines in the "ind" file, replacing population
    # labels with their substitutes
    for (merge_into in names(merge)) {
        regex <- paste0("(", paste(merge[[merge_into]], collapse="|"), ")$")
        lines <- str_replace(lines, regex, merge_into)
    }

    writeLines(lines, modified_file)
}



# Reading EIGENSTRAT files --------------------------------------------------


#' Read an EIGENSTRAT 'ind' file.
#'
#' @param file Path to the file.
#'
#' @return Data frame with the sample identifier, sex and label
#'     columns (columns defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
read_ind <- function(file) {
    read_table2(file, col_names=c("id", "sex", "label"))
}


#' Read an EIGENSTRAT 'snp' file.
#'
#' @param file Path to the file.
#'
#' @return Data frame with information about each SNP (columns defined
#'     by the EIGENSTRAT format).
#'
#' @export
#' @import readr
read_snp <- function(snp_file) {
    read_table2(snp_file, col_names=c("id", "chrom", "gen", "pos", "ref", "alt"))
}


#' Read an EIGENSTRAT 'geno' file.
#'
#' @param file Path to the geno file.
#' @param ind_file Path to the ind file to read sample names from.
#'
#' @return Data frame with columns containing "genotypes" of each
#'     sample (0/1/9 as defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
read_geno <- function(file, ind_file=NULL) {
    if (!is.null(ind_file)) {
        inds <- read_ind(ind_file)$id
    } else {
        inds <- NULL
    }

    # get the number of samples in the geno file
    n <- nchar(readLines(file, 1))
    read_fwf(file, col_positions=fwf_widths(rep(1, n), inds))
}


#' Read a tripplet of EIGENSTRAT (geno/snp/ind files) files.
#'
#' @param file Path to the file.
#'
#' @return List of three data frames (one element for geno/snp/ind).
#'
#' @export
read_eigenstrat <- function(prefix=NULL) {
    list(
        geno=read_geno(paste0(prefix, ".geno"), paste0(prefix, ".ind")),
        snp=read_snp(paste0(prefix, ".snp")),
        ind=read_ind(paste0(prefix, ".ind"))
    )
}



# Writing EIGENSTRAT files --------------------------------------------------


#' Write an EIGENSTRAT 'ind' file.
#'
#' @param file Path to the file.
#' @param ind Data frame with the sample identifier, sex and label
#'     columns (columns defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
write_ind <- function(ind_file, df) {
    write_tsv(df, ind_file, col_names=FALSE)
}


#' Write an EIGENSTRAT 'snp' file.
#'
#' @param file Path to the file.
#' @param snp Data frame with information about each SNP (columns
#'     defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
write_snp <- function(snp_file, df) {
    write_tsv(df, snp_file, col_names=FALSE)
}


#' Write an EIGENSTRAT 'geno' file.
#'
#' @param file Path to the file.
#' @param geno Data frame with columns containing "genotypes" of each
#'     sample (0/1/9 as defined by the EIGENSTRAT format).
#'
#' @export
#' @import readr
write_geno <- function(geno_file, df) {
    writeLines(apply(df, 1, paste, collapse=""), con=geno_file)
}


#' Write a tripplet of EIGENSTRAT (geno/snp/ind files) files.
#'
#' @param prefix Prefix of the geno/snp/ind files (can
#'     include the path).
#' @param ind data.frame with data in a 'ind' format
#' @param snp data.frame with data in a 'snp' format
#' @param geno data.frame with data in a 'geno' format
#'
#' @return List of three data frames (one element for geno/snp/ind).
#'
#' @export
write_eigenstrat <- function(prefix, ind, snp, geno) {
    write_ind(paste0(prefix, ".ind"), ind)
    write_snp(paste0(prefix, ".snp"), snp)
    write_geno(paste0(prefix, ".geno"), geno)
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
#' @param pos_file Path to the 2 column position file to intersect with.
#' @param complement Perform an intersect or a complement operation?
#'
#' @import readr dplyr
#' @export
subset_sites <- function(prefix, out_prefix, bed_file=NULL, pos_file=NULL, complement=FALSE) {
    if (!is.null(bed_file)) {
        coords <- read_table2(bed_file, col_names=c("chrom", "start", "end")) %>% select(chrom, end)
    } else if(!is.null(pos_file)) {
        coords <- read_table2(pos_file, col_names=c("chrom", "pos"))
    }

    geno <- read_geno(paste0(prefix, ".geno"))
    snp <- read_snp(paste0(prefix, ".snp"))
    combined <- bind_cols(snp, geno)

    # determine which function to call on the coordinates
    fun <- ifelse(complement, anti_join, inner_join)
    combined_subset <- fun(combined, coords)

    # write the new snp file
    write_tsv(select(combined_subset, id:ref), path=paste0(out_prefix, ".snp"), col_names=FALSE)
    # write the new geno file
    writeLines(apply(select(combined_subset, -(id:ref)), 1, paste, collapse=""),
               con=paste0(out_prefix, ".geno"))
    # write the new ind file
    file.copy(from=paste0(prefix, ".ind"), to=paste0(out_prefix, ".ind"))
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
#' @import readr
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
#' @import readr
#' @export
vcf_to_eigenstrat <- function(vcf_file, prefix) {
    vcf <- read_tsv(vcf_file, comment="##") %>%
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
