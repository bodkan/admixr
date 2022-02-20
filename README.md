

# admixr

<!-- badges: start -->
[![R-CMD-check](https://github.com/bodkan/admixr/workflows/R-CMD-check/badge.svg)](https://github.com/bodkan/admixr/actions)
[![Coverage status](https://codecov.io/gh/bodkan/admixr/branch/main/graph/badge.svg)](https://codecov.io/github/bodkan/admixr?branch=main)
<!-- badges: end -->

## Introduction

[ADMIXTOOLS](https://github.com/DReichLab/AdmixTools/) is a widely used
software package for calculating admixture statistics and testing population
admixture hypotheses.

A typical ADMIXTOOLS workflow often involves a combination of `sed`/`awk`/shell
scripting and manual editing to create different configuration files. These are
then passed as command-line arguments to one of ADMIXTOOLS commands, and
control how to run a particular analysis. The results are then redirected to
another file, which has to be parsed by the user to extract values of interest,
often using command-line utilities again or (worse) by manual copy-pasting.
Finally, the processed results are analysed in R, Excel or another program.

This workflow can be a little cumbersome, especially if one wants to explore many
hypotheses involving different combinations of populations. Most importantly,
however, it makes it difficult to follow the rules of best practice for
reproducible science, as it is nearly impossible to construct reproducible
automated "pipelines".

This R package makes it possible to perform all stages of an ADMIXTOOLS
analysis entirely from R. It provides a set of convenient functions that
completely remove the need for "low level" configuration of individual
ADMIXTOOLS programs, allowing users to focus on the analysis itself.

## How to cite

_admixr_ is now published
as an [_Application Note_](https://doi.org/10.1093/bioinformatics/btz030) in the journal Bioinformatics. If you use it in your work, please cite the paper!

## Installation instructions

##### Latest stable version

The package is available [on
CRAN](https://cran.r-project.org/package=admixr). You can install it
simply by running


```r
install.packages("admixr")
```

from your R session. This the recommended procedure for most users.

##### Development version

To install the development version from Github (which might be
slightly ahead in terms of new features and bugfixes compared to the
stable release on CRAN), you need the package `devtools`. You can run:


```r
install.packages("devtools")
devtools::install_github("bodkan/admixr")
```

##### Installing ADMIXTOOLS

In order to use the `admixr` package, you need a working installation
of ADMIXTOOLS. You can find installation instructions
[here](https://github.com/DReichLab/AdmixTools/blob/master/README.INSTALL).

Furthermore, you also need to make sure that R can find ADMIXTOOLS
binaries on the `$PATH`. You can achieve this by specifying
`PATH=<path to the location of ADMIXTOOLS programs>` in the
`.Renviron` file.

## Example

This is all the code that you need to perform ADMIXTOOLS analyses using this
package! No shell scripting, no copy-pasting and manual editing of text files.
The only thing you need is a working ADMIXTOOLS installation and a path to
EIGENSTRAT data (a trio of ind/snp/geno files), which we call `prefix` here.


```r
library(admixr)

# download a small testing dataset to a temporary directory and
# process it for use in R
snp_data <- eigenstrat(download_data())

result <- d(
  W = c("French", "Sardinian"), X = "Yoruba", Y = "Vindija", Z = "Chimp",
  data = snp_data
)

result
#> # A tibble: 2 × 10
#>   W         X      Y       Z          D  stderr Zscore  BABA  ABBA  nsnps
#>   <chr>     <chr>  <chr>   <chr>  <dbl>   <dbl>  <dbl> <dbl> <dbl>  <dbl>
#> 1 French    Yoruba Vindija Chimp 0.0313 0.00693   4.51 15802 14844 487753
#> 2 Sardinian Yoruba Vindija Chimp 0.0287 0.00679   4.22 15729 14852 487646
```

Note that a single call to the `d` function generates all required intermediate
config and population files, runs ADMIXTOOLS, parses its log output and returns
the result as a `data.frame` object. It does all of this behind the scenes,
without the user having to deal with low-level technical details.

## More information

To see many more examples of admixr in action, please check out the
[tutorial vignette](https://bodkan.net/admixr/articles/tutorial.html).

If you want to stay updated on new _admixr_ development, follow me on
[Twitter](https://www.twitter.com/fleventy5).
