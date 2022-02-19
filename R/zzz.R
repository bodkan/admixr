.onAttach <- function(libname, pkgname) {
  ## check for presence of an ADMIXTOOLS command on user's PATH and display
  ## a warning if it's not present
  cmds  <- c("qp3Pop", "qpDstat", "qpF4ratio", "qpAdm", "qpWave", "qpfstats", "mergeit")
  path_check <- all(Sys.which(cmds) != "")
  if (!path_check) {
    packageStartupMessage(
      "ADMIXTOOLS programs required by admixr are not in your $PATH variable.\n",
      "At the very least, the following programs should be available:\n    ",
      paste(cmds, collapse = ", "),
      "\nOne way to solve this is to modify the $PATH variable in your ~/.Renviron\n",
      "file so that it points to the directory containing the ADMIXTOOLS programs.\n\n",
      "Current contents of the $PATH:\n", Sys.getenv("PATH", "\n")
    )
  }
}
