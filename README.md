

# _admixr_&mdash;interactive R interface for ADMIXTOOLS

[![CRAN-version](https://www.r-pkg.org/badges/version/admixr)](https://cran.r-project.org/package=admixr) [![CRAN-downloads](https://cranlogs.r-pkg.org/badges/grand-total/admixr)](https://cran.r-project.org/package=admixr) [![R-CMD-check](https://github.com/bodkan/admixr/workflows/R-CMD-check/badge.svg)](https://github.com/bodkan/admixr/actions) [![Binder](http://mybinder.org/badge.svg)](http://beta.mybinder.org/v2/gh/bodkan/admixr/main?urlpath=rstudio) [![Coverage status](https://codecov.io/gh/bodkan/admixr/branch/main/graph/badge.svg)](https://app.codecov.io/github/bodkan/admixr?branch=main)

<!-- badges: end -->

## What is _admixr_?

The _admixr_ package provides a convenient R interface to
[ADMIXTOOLS](https://github.com/DReichLab/AdmixTools/), a widely used
software package for calculating admixture statistics and testing population
admixture hypotheses.

A typical ADMIXTOOLS workflow often involves a combination of `sed`/`awk`/shell
scripting and manual editing to create different configuration files. These are
then passed as command-line arguments to one of ADMIXTOOLS commands, and
control how to run a particular analysis. The results of such computation are
then usually redirected to another file, which needs to be parsed by the user
to extract values of interest, often using command-line utilities again or by
manual copy-pasting, and finally analysed in R, Excel or another program.

This workflow can be a little cumbersome, especially if one wants to explore many
hypotheses involving different combinations of populations or data filtering
strategies. Most importantly, it makes it difficult to follow the rules of best
practice for reproducible science, especially given the need for manual
intervention on the command-line or custom shell scripting to orchestrate more
complex pipelines.

_admixr_ makes it possible to perform all stages of an ADMIXTOOLS analysis entirely
from R. It provides a [set of convenient functions](https://bodkan.net/admixr/reference/index.html)
that completely remove the need for "low-level" configuration of individual ADMIXTOOLS
programs, allowing users to focus on the analysis itself.

## How to cite

_admixr_ is now published
as an [_Application Note_](https://doi.org/10.1093/bioinformatics/btz030) in the
journal Bioinformatics. If you use it in your work, please cite the paper! You
will join an [excellent company](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=13286994334855947290)
of papers who have used it to do amazing research. 🙂

## Installation instructions

#### Browser-based RStudio session

You can try out _admixr_ without installation directly in your browser! Simply
click on [![Binder](http://mybinder.org/badge.svg)](http://beta.mybinder.org/v2/gh/bodkan/admixr/main?urlpath=rstudio)
and after a short moment you will get a Binder RStudio could session running
in your web browser. However, please note that Binder's computational resources
are extremely limited so you might run into issues if you try to run extremely
resource-intensive computations.

#### Latest stable version

The package is available [on
CRAN](https://cran.r-project.org/package=admixr). You can install it
simply by running


``` r
install.packages("admixr")
```

from your R session. This the recommended procedure for most users.

#### Development version

To install the development version from Github (which might be
slightly ahead in terms of new features and bugfixes compared to the
stable release on CRAN), you need the R package _devtools_. You can run:


``` r
install.packages("devtools")
devtools::install_github("bodkan/admixr")
```

#### Installing ADMIXTOOLS

In order to use the _admixr_ package, you need a working installation
of ADMIXTOOLS. You can find installation instructions
[here](https://github.com/DReichLab/AdmixTools/blob/master/README.INSTALL).

Furthermore, you also need to make sure that R can find ADMIXTOOLS
binaries on the `$PATH`. You can achieve this by specifying
`PATH=<path to the location of ADMIXTOOLS programs>` in the
`.Renviron` file in your home directory. If R cannot find ADMIXTOOLS utilities,
you will get a warning upon loading `library(admixr)` in your R session.

**In terms of the version of ADMIXTOOLS required by _admixr_, there is no
formal requirement. However, please note that the package is generally developed
and tested against the latest ADMIXTOOLS release (which is currently version 8.0.2).**

## Example analysis

This is all the code that you need to perform ADMIXTOOLS analyses using this
package! No shell scripting, no copy-pasting and manual editing of text files.
The only thing you need is a working ADMIXTOOLS installation and a path to
EIGENSTRAT data (a trio of ind/snp/geno files), which we call `prefix` here.


``` r
library(admixr)

# download a small testing dataset to a temporary directory and process it for use in R
snp_data <- eigenstrat(download_data())

result <- d(
  W = c("French", "Sardinian"), X = "Yoruba", Y = "Vindija", Z = "Chimp",
  data = snp_data
)

result
```

Note that a single call to the `d` function generates all required intermediate
config and population files, runs ADMIXTOOLS, parses its log output and returns
the result as a `data.frame` object with the D statistics results. It does all of
this behind the scenes, without the user having to deal with low-level technical
details.

## More information

To see many more examples of admixr in action, please check out the
[tutorial vignette](https://bodkan.net/admixr/reference/index.html).

## Is _admixr_ related to ADMIXTOOLS 2?

Recently, a new R package called [ADMIXTOOLS 2](https://uqrmaie1.github.io/admixtools/)
appeared on the horizon, offering a re-implementation of several features of the
original ADMIXTOOLS suite of command-line programs.

The _admixr_ project is not related to that initiative. It is not a pre-cursor to it, nor
it is superseeded by it. I have never used ADMIXTOOLS 2 myself, but from the looks of it
it seems to offer some very interesting features for fitting complex admixture graphs,
which is certainly something which _admixr_ does not do.

**The bottom-line is this:** as long as the [original ADMIXTOOLS](https://github.com/DReichLab/AdmixTools)
continues to be developed and maintained, _admixr_ remains relevant and useful and
will continue to be supported. ADMIXTOOLS is one of the most battle-tested pieces
of software in population genetics&mdash;if you're happy with the set of features
it provides and if you're happy with _admixr_ itself, there is no real reason
to move away from either of them.
