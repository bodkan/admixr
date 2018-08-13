
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build
Status](https://travis-ci.org/bodkan/admixr.svg?branch=master)](https://travis-ci.org/bodkan/admixr)
[![Coverage
status](https://codecov.io/gh/bodkan/admixr/branch/master/graph/badge.svg)](https://codecov.io/github/bodkan/admixr?branch=master)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

# admixr

***“Make ADMIXTOOLS fun again\!”***

**TL;DR** Perform ADMIXTOOLS analyses directly from R, without ever
having to deal with par/pop/log files again.

## Introduction

[ADMIXTOOLS](http://www.genetics.org/content/192/3/1065) is a widely
used software package for calculating statistics and testing hypotheses
about population admixture. However, although very powerful and
comprehensive, it leaves a lot to be desired in terms of user
experience. A typical ADMIXTOOLS workflow usually involves a combination
of `sed`/`awk`/shell scripting and manual editing to create different
configuration files. These files are then supplied by the user as
command-line arguments to one of the ADMIXTOOLS commands, and describe
how to run a particular analysis. The results of each analysis are then
usually redirected to another file. Afterwards, the user needs to
extract values of interest from this file (which can be somewhat complex
and full of redundant information), most likely using more command-line
scripting or manual copy-pasting, and analyse them in R or Excel. This
workflow is quite cumbersome, especially if one wants to test many
hypotheses involving different combinations of populations. Most
importantly, however, it makes it difficult to follow good practices of
reproducible science.

This R package makes it possible to perform all stages of ADMIXTOOLS
analysis entirely from within R. It provides a set of convenient
functions that completely abstract away the need for “low level”
configuration of individual ADMIXTOOLS programs, allowing the user to
focus on the analysis itself.

## Installation instructions

To install `admixr` from Github you need to install the package
`devtools` first. You can simply run:

    install.packages("devtools")
    devtools::install_github("bodkan/admixr")

## Example

This is all the code that you need to perform ADMIXTOOLS analyses using
this package\! No shell scripting, no copy-pasting and manual editing of
text files. The only thing you need is a working ADMIXTOOLS installation
and a path to EIGENSTRAT data (a trio of ind/snp/geno files).

``` r
library(admixr)

res <- d(
  W = "Yoruba", X = "French", Y = c("Han", "Japanese"), Z = "Uygur",
  prefix = "/Users/martin_petr/local/AdmixTools-5.0/data/allmap"
)

res
#>        W      X        Y     Z      D   stderr Zscore  BABA  ABBA  nsnps
#> 1 Yoruba French      Han Uygur 0.0410 0.001271 32.300 30753 28328 621026
#> 2 Yoruba French Japanese Uygur 0.0402 0.001295 31.016 30690 28320 621026
```

Note that a single call to the `d` function generates all required
intermediate config and population files, runs ADMIXTOOLS, parses its
log output and returns results as a `data.frame` object. It does so all
behind the scenes without the user having to deal with “low-level”
technical details.
