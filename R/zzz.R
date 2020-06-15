.onAttach <- function(libname, pkgname) {
  # check for presence of an ADMIXTOOLS command on user's PATH and display
  # a warning if it's not present
  path_check <- system("command -v qpDstat", ignore.stdout = TRUE)
  if (path_check != 0) {
    packageStartupMessage(
      "ADMIXTOOLS binaries could not be found in your $PATH.\n",
      "Please make sure you have ADMIXTOOLS compiled and that all commands are in $PATH.",
      "\nThen make sure to modify the $PATH variable in your .Renviron file",
      " so that R can find the ADMIXTOOLS binaries (you can run 'echo \"PATH=$PATH\" >> .Renviron' in the shell)."
    )
  }
}
