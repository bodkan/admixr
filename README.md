# admixr

_**Make ADMIXTOOLS fun again!**_

This package makes it possible to perform all the stages of ADMIXTOOLS
analysis without leaving R. A typical ADMIXTOOLS workflow usually
involves using a combination of `sed`, `awk` and manual editing to
create the parameter and population configuration files that specify
which ADMIXTOOLS command to run and how to run it. Afterwards, the
user needs to extract the results from an output file (which can be
somewhat complex and full of redundant information), most likely using
command-line utilities again, and analyse them in R or Excel. This can
be cumbersome especially if one wants to run many analyses which test
different hypotheses that involve various combinations of
populations. This package provides a set of functions that make it
easy to automatically generate the required configuration files, as
well as a set of wrapper functions that run different ADMIXTOOLS
commands based on these automatically generated configuration
files. Finally, it provides a simple way to parse ADMIXTOOLS output
files, extracting only the useful information and returning a simple
dataframe of all results for immediate plotting and statistical
analysis.

## Installation instructions

To install `admixr` from Github you need to install the package `devtools` first. After you have done that, it is simply a matter of:

```
devtools::install_github("bodkan/admixr")
library(admixr)
```

## Warning!

This package is not even in an alpha stage yet and the API is still changing. It only supports `qpDstat` and `qpF4ratio` at this moment, although all ADMIXTOOLS commands will eventually be supported.

You probably shouldn't use it yet. :)
