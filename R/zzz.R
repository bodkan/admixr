.onAttach <- function(libname, pkgname) {
  ## check for presence of an ADMIXTOOLS command on user's PATH and display
    ## a warning if it's not present
    cmds  <- c("qp3Pop", "qpDstat", "qpF4ratio", "qpAdm", "qpWave")
    path_check <- all(Sys.which(cmds) != "")
    if (!path_check) {
        packageStartupMessage(
            "Not all ADMIXTOOLS binaries could be found in your $PATH variable.\n",
            "At the very least, the following programs are recommended:\n    ",
            paste(cmds, collapse = ", "),
            "\nMake sure to modify the $PATH variable in your .Renviron file so \n",
            "that it points to the directory containing the ADMIXTOOLS programs.\n")
    }
}
