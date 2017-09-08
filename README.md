# admixr - A set of convenient functions for working with ADMIXTOOLS

This package makes it possible to perform all the stages of ADMIXTOOLS
analysis without leaving R. A typical ADMIXTOOLS workflow usually
involves using a combination of sed, awk and manual editing to create
a set of so called parameter and population files, which can be
cumbersome especially in case of running many analyses, involving
different combinations of populations. This package provides a
convenient set of functions to use the scripting capabilities of R to
automatically generate the required configuration files from R, as
well as a set of wrapper function to use those configuration files to
run ADMIXTOOLS commands. Finally, it makes it possible to parse the
output files of ADMIXTOOLS programs (which can be somewhat comples and
full of redundant information), returning to the user a simple
dataframe of all results, which can be immediately used for plotting
and running statistical analyses.