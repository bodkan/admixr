#' Generate a parameter file.
#'
#' @param files List of filenames of the population file, parameter
#'     file and log file.
#' @param prefix Prefix of the geno/snp/ind files (can include the
#'     path). If specified, geno_file/snp_file/ind_file will be
#'     ignored.
#' @param geno_file Path to the genotype file.
#' @param snp_file Path to the snp file.
#' @param ind_file Path to the ind file.
#' @param badsnp_file SNP file with information about ignored sites.
#' @export
create_par_file <- function(files,
                            prefix=NULL,
                            geno_file=NULL, snp_file=NULL, ind_file= NULL,
                            badsnp_file=NULL) {
    if (all(is.null(c(prefix, geno_file, snp_file, ind_file)))) {
        stop("Prefix of EIGENSTRAT files or the paths to individual geno/snp/ind files must be specified")
    }

    # if the user specified EIGENSTRAT prefix, set only paths to unspecified geno/snp/ind files
    if (!is.null(prefix)) {
        if (is.null(geno_file)) geno_file <- paste0(prefix, ".geno")
        if (is.null(snp_file)) snp_file <- paste0(prefix, ".snp")
        if (is.null(ind_file)) ind_file <- paste0(prefix, ".ind")
    }

    writeLines(sprintf("genotypename: %s\nsnpname: %s\nindivname: %s\n",
                       geno_file, snp_file, ind_file),
               con=files$par_file)

    if (!is.null(files[["pop_file"]])) {
        write(sprintf("popfilename: %s\n", files$pop_file), file=files$par_file, append=TRUE)
    } else if (!is.null(files$popleft) & !is.null(files$popright)) {
        write(sprintf("popleft: %s", files$popleft), file=files$par_file, append=TRUE)
        write(sprintf("popright: %s", files$popright), file=files$par_file, append=TRUE)
    }

    if (!is.null(badsnp_file)) {
        write(sprintf("badsnpname: %s", badsnp_file), file=files$par_file, append=TRUE)
    }
}


# Generate a file with populations for a qpF4ratio run.
create_qpF4ratio_pop_file <- function(X, A, B, C, O, file) {
    lines <- sprintf("%s %s : %s %s :: %s %s : %s %s", A, O, X, C, A, O, B, C)
    writeLines(lines, file)
}


# Generate a file with populations for a qpDstat run.
create_qpDstat_pop_file <- function(W=NULL, X=NULL, Y=NULL, Z=NULL, file) {
    lines <- c()
    for (w in W) for (x in X) for (y in Y) for (z in Z) {
        lines <- c(lines, sprintf("%s %s %s %s", w, x, y, z))
    }
    writeLines(lines, file)
}


# Generate a file with populations for a qp3Pop run.
create_qp3Pop_pop_file <- function(A, B, C, file) {
    lines <- c()
    for (a in A) for (b in B) for (c in C) {
        lines <- c(lines, sprintf("%s %s %s", a, b, c))
    }
    writeLines(lines, file)
}


# Generate a file with populations for a qpAdm run.
#
# L, R - Sets of left (U) and right (R) populations using the
#   terminology of Haak et al., 2012 (Supplementary Information 10 on
#   page 128).
# files - A list that must contain "popleft" and "popright" elements,
#   which describe the paths to files containing "left" and "right"
#   populations (one population per line).
create_qpAdm_pop_files <- function(L, R, files) {
    writeLines(L, con=files[["popleft"]])
    writeLines(R, con=files[["popright"]])
}