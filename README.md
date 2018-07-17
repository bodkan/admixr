
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build
Status](https://travis-ci.org/bodkan/admixr.svg?branch=master)](https://travis-ci.org/bodkan/admixr)

# admixr

***“Make ADMIXTOOLS fun again\!”***

**TL;DR** Perform ADMIXTOOLS analyses directly from R, without having to
touch ADMIXTOOLS/parfiles/popfiles/logfiles ever again.

## Introduction

[ADMIXTOOLS](http://www.genetics.org/content/192/3/1065) is a popular
and widely used software package for calculating various population
genetic statistics and for testing population genetics hypotheses.
However, although very powerful and comprehensive, it leaves a lot to be
desired in terms of user experience.

This R package aims to make it possible to perform all stages of
ADMIXTOOLS analysis entirely from within R. A typical ADMIXTOOLS
workflow usually involves using a combination of `sed`/`awk`/shell
scripting and manual editing to create parameter and population
configuration files. These files are then supplied as arguments to an
ADMIXTOOLS command and describe how to run a particular analysis, with
results saved to another file. Afterwards, the user needs to extract
values of interest from this file (which can be somewhat complex and
full of redundant information), most likely using more command-line
scripting, and analyse them in R or Excel. This can be cumbersome
especially if one wants to run a large number of analyses that can
involve many combinations of populations.

This package provides a set of functions that completely abstract away
the “low level” configuration of ADMIXTOOLS programs, making it possible
to focus on the analysis itself. It achieves this by automatically
generating all configuration files, running the commands and parsing
their outputs behind the scenes, extracting only useful information and
always returning a simple dataframe for downstream analyses.

## Installation instructions

To install `admixr` from Github you need to install the package
`devtools` first. Then it’s a simple matter of running:

    install.packages("devtools")
    devtools::install_github("bodkan/admixr")

## Warning\!

**This package is an alpha stage software and its API is still
changing\!** It only supports only a subset of available ADMIXTOOLS
commands, although all ADMIXTOOLS commands will eventually be supported.

## Example

The following assumes that we have a working installation of ADMIXTOOLS
in `~/local/AdmixTools-5.0/`. This example requires a small EIGENSTRAT
data set downloaded for running testing script of the ADMIXTOOLS
package.

Note that a single call to our `qpDstat` function generates all required
intermediate config and population files, runs ADMIXTOOLS, parses its
log output and returns results as a `data.frame` object. It does so all
behind the scenes without the user having to deal with “low-level”
technical details.

``` r
library(admixr)

res <- qpDstat(
  W = "Yoruba", X = "French", Y = c("Han", "Japanese"), Z = "Uygur",
  prefix = "/Users/martin_petr/local/AdmixTools-5.0/data/allmap"
)

res
#> # A tibble: 2 x 10
#>   W      X      Y        Z          D  stderr Zscore  BABA  ABBA n_snps
#>   <chr>  <chr>  <chr>    <chr>  <dbl>   <dbl>  <dbl> <int> <int>  <int>
#> 1 Yoruba French Han      Uygur 0.041  0.00127   32.3 30753 28328 621026
#> 2 Yoruba French Japanese Uygur 0.0402 0.00130   31.0 30690 28320 621026
```
