# Generate a parameter file.
create_par_file <- function(files, data, params) {
  geno <- data$geno
  snp <- data$snp
  ind <- ifelse(is.null(data$group), data$ind, data$group)
  exclude <- data$exclude

  sprintf("genotypename: %s\nsnpname: %s\nindivname: %s\n", geno, snp, ind) %>%
    writeLines(con = files$par_file)

  if (!is.null(files[["pop_file"]])) {
    write(sprintf("popfilename: %s\n", files$pop_file), file = files$par_file, append = TRUE)
  } else if (!is.null(files$popleft) & !is.null(files$popright)) {
    write(sprintf("popleft: %s", files$popleft), file = files$par_file, append = TRUE)
    write(sprintf("popright: %s", files$popright), file = files$par_file, append = TRUE)
  }

  if (!is.null(exclude)) {
    write(sprintf("badsnpname: %s", exclude), file = files$par_file, append = TRUE)
  }

  # write user-specified parameters
  for (par in names(params)) {
    write(sprintf("%s: %s\n", par, params[[par]]), file = files$par_file, append = TRUE)
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


# Generate a file with populations for a qpDstat run based on given
# list of quartets.
create_qpDstat_pop_file_quartets <- function(quartets, file) {
  lines <- c()
  for (q in quartets) {
    lines <- c(lines, sprintf("%s %s %s %s", q[1], q[2], q[3], q[4]))
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
create_leftright_pop_files <- function(left, right, files) {
  writeLines(left, files$popleft)
  writeLines(right, files$popright)
}

