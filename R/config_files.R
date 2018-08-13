# Generate a parameter file.
create_par_file <- function(files,
                            prefix = NULL,
                            geno_file = NULL, snp_file = NULL, ind_file = NULL,
                            badsnp_file = NULL) {
    if (all(is.null(c(prefix, geno_file, snp_file, ind_file)))) {
        stop("EIGENSTRAT prefix geno/snp/ind paths must be specified")
    }

    # if the user specified EIGENSTRAT prefix, set only paths to unspecified
    # geno/snp/ind files
    if (!is.null(prefix)) {
      prefix <- path.expand(prefix)
      if (is.null(geno_file)) geno_file <- paste0(prefix, ".geno")
      if (is.null(snp_file)) snp_file <- paste0(prefix, ".snp")
      if (is.null(ind_file)) ind_file <- paste0(prefix, ".ind")
    } else {
      ind_file <- path.expand(ind_file)
      snp_file <- path.expand(snp_file)
      geno_file <- path.expand(geno_file)
    }

    writeLines(sprintf("genotypename: %s\nsnpname: %s\nindivname: %s\n",
                       geno_file, snp_file, ind_file),
               con = files$par_file)

    if (!is.null(files[["pop_file"]])) {
        write(sprintf("popfilename: %s\n", files$pop_file), file = files$par_file, append = TRUE)
    } else if (!is.null(files$popleft) & !is.null(files$popright)) {
        write(sprintf("popleft: %s", files$popleft), file = files$par_file, append = TRUE)
        write(sprintf("popright: %s", files$popright), file = files$par_file, append = TRUE)
    }

    if (!is.null(badsnp_file)) {
        write(sprintf("badsnpname: %s", badsnp_file), file = files$par_file, append = TRUE)
    }
}


# Generate a file with populations for a qpF4ratio run.
create_qpF4ratio_pop_file <- function(X, A, B, C, O, file) {
    lines <- sprintf("%s %s : %s %s :: %s %s : %s %s", A, O, X, C, A, O, B, C)
    writeLines(lines, file)
}


# Generate a file with populations for a qpDstat run.
create_qpDstat_pop_file <- function(W = NULL, X = NULL, Y = NULL, Z = NULL, file) {
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

# Generate a file with populations for a qpDstat run.
create_qpAdm_pop_files <- function(left, right, files) {
  writeLines(left, files$popleft)
  writeLines(right, files$popright)
}