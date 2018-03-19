# admixr

[![Travis-CI Build Status](https://travis-ci.org/bodkan/admixr.svg?branch=master)](https://travis-ci.org/bodkan/admixr)

_**"Make ADMIXTOOLS fun again!"**_

[ADMIXTOOLS](http://www.genetics.org/content/192/3/1065) is a popular
and widely used software package for calculating various population
genetic statistics and hypotheses testing. However, although incredibly
powerful and comprehensive, in the user experience department it leaves
a lot to be desired.

This R package aims to allow to perform all stages of ADMIXTOOLS
analysis without leaving R. A typical ADMIXTOOLS workflow usually
involves using a combination of `sed`, `awk` and manual editing to
create parameter and population configuration files for a given ADMIXTOOLS
command to run and describing how to run it. Afterwards, the
user needs to extract results from an output file (which can be
somewhat complex and full of redundant information), most likely using
command-line utilities, and analyse them in R or Excel. This can
be cumbersome especially if one wants to run many analyses which test
different hypotheses that involve many combinations of populations.

This package provides a set of functions that make it easy to automatically
generate  required configuration files, as well as a set of wrapper
functions that run different ADMIXTOOLS commands based on these
files. Finally, it provides a consistent way to parse ADMIXTOOLS output
files, extracting only the useful information and returning a simple
dataframe of all results for downstream analyses and plotting.

## Installation instructions

To install `admixr` from Github you need to install the package `devtools`
first. After you have done that, it is simply a matter of:

```
devtools::install_github("bodkan/admixr")
library(admixr)
```

## Warning!

This package is not even in an alpha stage yet and the API is still changing.
It only supports `qpDstat` and `qpF4ratio` at this moment, although all ADMIXTOOLS
commands will eventually be supported.

**You probably shouldn't use it yet. :)**
