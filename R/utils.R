#' Generate a file with populations for a qpF4ratio run.
#'
#' @description For details about the format of the population list
#'     file see documentation of the ADMIXTOOLS package.
#' 
#' @param X, A, B, C, O Population names, using the terminology of
#'     Patterson et al., 2012
#' @param file Path to the poplist file that will be generated.
#' @export
create_qpF4ratio_pop_file <- function(X, A, B, C, O, file) {
    lines <- sprintf("%s %s : %s %s :: %s %s : %s %s", A, O, X, C, A, O, B, C)
    writeLines(lines, file)
}


#' Generate a file with populations for a qpDstat run.
#'
#' @description For details about the format of the population list
#'     file see documentation of the ADMIXTOOLS package.
#'
#' @param W, X, Y, Z Population names, using the terminology of
#'     Patterson et al., 2012. It is possible to specify vectors of
#'     population identifiers as W, X, etc. That way, all possible
#'     combinations of specified W, X, Y and Z quadruples will be
#'     generated.
#' @param poplist List of populations. It will be saved as a file, one
#'     population per line.
#' @param file Path to the poplist file that will be generated.
#' @export
create_qpDstat_pop_file <- function(W=NULL, X=NULL, Y=NULL, Z=NULL, file) {
    lines <- c()
    for (w in W) for (x in X) for (y in Y) for (z in Z) {
        lines <- c(lines, sprintf("%s %s %s %s", w, x, y, z))
    }
    writeLines(lines, file)
}


#' Generate a file with populations for a qp3Pop run.
#'
#' @description For details about the format of the population list
#'     file see documentation of the ADMIXTOOLS package.
#'
#' @param A, B, C Population names, using the terminology of Patterson
#'     et al., 2012. It is possible to specify vectors of population
#'     identifiers as A, B, C. That way, all possible triplets of
#'     populations will be generated.
#' @param file Path to the poplist file that will be generated.
#' @export
create_qp3Pop_pop_file <- function(A, B, C, file) {
    lines <- c()
    for (a in A) for (b in B) for (c in C) {
        lines <- c(lines, sprintf("%s %s %s", a, b, c))
    }
    writeLines(lines, file)
}


#' Generate a file with populations for a qpAdm run.
#'
#' @description For details about the format of the population list
#'     file see documentation of the ADMIXTOOLS package.
#'
#' @param L, R Sets of left (U) and right (R) populations using the
#'     terminology of Haak et al., 2012 (Supplementary Information 10
#'     on page 128).
#' @param files A list that must contain "popleft" and "popright" elements, which describe the paths to files containing "left" and "right" populations (one population per line).
#'
#' @export
create_qpAdm_pop_files <- function(L, R, files) {
    writeLines(L, con=files[["popleft"]])
    writeLines(R, con=files[["popright"]])
}


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


#' Run a given ADMIXTOOLS command.
#'
#' @param tool Which ADMIXTOOLS command to run. At the moment, only
#'     qpDstat and qpF4ratio are supported.
#' @param par_file Path to the parameter file.
#' @param log_file Path to the output file. If NULL, output will be
#'     printed to stdout.
#' @export
run_cmd <- function(cmd, par_file, log_file) {
    if (!cmd %in% c("qpDstat", "qpF4ratio", "qp3Pop", "qpAdm")) {
        stop("ADMIXTOOLS command '", cmd, "' is not supported or does not exist")
    }

    system(paste(cmd, "-p", par_file, ">", log_file))
}


# Create either specified or a temporary directory.
get_dir <- function(dir_name=NULL) {
    if (!is.null(dir_name)) {
        dir.create(dir_name, showWarnings=FALSE)
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
        pop_file=file.path(directory, paste0(prefix, ".pop")),
        par_file=file.path(directory, paste0(prefix, ".par")),
        log_file=file.path(directory, paste0(prefix, ".log"))
    )
}

# Check for the presence of a given set of labels in an 'ind' file.
# Fail if there a sample was not found.
check_presence <- function(labels, prefix=NULL, ind=NULL) {
    if (!is.null(ind)) {
        path <- ind
    } else {
        path <- paste0(prefix, ".ind")
    }

    not_present <- setdiff(labels, suppressMessages(read_ind(path)$label))
    if (length(not_present) > 0) {
        stop("The following samples are not present in '", ind, "': ",
             paste(not_present, collapse=", "))
    }
}


# Convert VCF-like GT string(s) into EIGENSTRAT genotypes.
# gt_to_eigenstrat(c(".|.", "./.", ".", "0|0", "0/0", "0", "0|1", "1|0", "0/1", "1|1", "1/1", "1"))
gt_to_eigenstrat <- function(gts) {
    eigen_gts <- gts %>%
        stringr::str_replace("^0$",   "2") %>%
        stringr::str_replace("^1$",   "0") %>%
        stringr::str_replace("^\\.\\|\\.$", "9") %>%
        stringr::str_replace("^\\./\\.$", "9") %>%
        stringr::str_replace("^\\.$",   "9") %>%
        stringr::str_replace("^0\\|0$", "2") %>%
        stringr::str_replace("^0/0$", "2") %>%
        stringr::str_replace("^0\\|1$", "1") %>%
        stringr::str_replace("^1\\|0$", "1") %>%
        stringr::str_replace("^0/1$", "1") %>%
        stringr::str_replace("^1\\|1$", "0") %>%
        stringr::str_replace("^1/1$", "0")

    eigen_gts
}

# Convert VCF-like GT string(s) into EIGENSTRAT genotypes.
# eigenstrat_to_gt(c(0, 1, 2, 9, 9, 2, 1, 0))
eigenstrat_to_gt <- function(eigenstrat_gts) {
    vcf_gts <- sapply(eigenstrat_gts, function(i) {
        if      (i == 0) "1/1"
        else if (i == 1) "0/1"
        else if (i == 2) "0/0"
        else             "./."
    })

    vcf_gts
}
