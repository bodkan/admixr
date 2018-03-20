[![Travis-CI Build Status](https://travis-ci.org/bodkan/admixr.svg?branch=master)](https://travis-ci.org/bodkan/admixr)

# admixr

_**"Make ADMIXTOOLS fun again!"**_

[ADMIXTOOLS](http://www.genetics.org/content/192/3/1065) is a popular
and widely used software package for calculating various population
genetic statistics and for testing population genetics hypotheses.
However, although very powerful and comprehensive, it leaves a lot to
be desired in terms of user experience.

This R package aims to make it possible to perform all stages of ADMIXTOOLS
analysis entirely from within R. A typical ADMIXTOOLS workflow usually
involves using a combination of `sed`/`awk`/shell scripting and manual
editing to create parameter and population configuration files. These files
are then supplied as arguments to an ADMIXTOOLS command and describe how to
run  a particular analysis, with results saved to another file. Afterwards,
the user needs to extract values of interest from this file (which can be
somewhat complex and full of redundant information), most likely using
more command-line scripting, and analyse them in R or Excel. This
can be cumbersome especially if one wants to run a large number of analyses
that can involve many combinations of populations.

This package provides a set of functions that completely abstract away the
"low level" configuration of ADMIXTOOLS programs, making it possible to
focus on the analysis itself. It achieves this by automatically
generating all configuration files, running the commands and parsing their
outputs behind the scenes, extracting only useful information and always
returning a simple dataframe for downstream analyses.

## Installation instructions

To install `admixr` from Github you need to install the package `devtools`
first. After you have done that, it is simply a matter of:

```
devtools::install_github("bodkan/admixr")
library(admixr)
```

## Warning!

**This package is an alpha stage software and its API is still changing!**
It only supports only a subset of available ADMIXTOOLS commands, although all
ADMIXTOOLS commands will eventually be supported.
