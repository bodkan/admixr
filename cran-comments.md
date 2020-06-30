## Resubmission

This is a third (re-)submission of the package to CRAN.

There were two new points raised in the last review:

1. After reading the first comment ("your examples in the
   documentation are all wrapped in 'do not run'"), I should clarify
   the following:

   The purpose of this package is to wrap around command-line programs
   in the ADMIXTOOLS suite for population genetics (a large C
   codebase, only available for unix/Linux). This package allows the
   user to do all the work in R (data processing/statistics/plotting)
   without ever touching bash/awk/sed/perl.

   All examples require working ADMIXTOOLS installation. This is why
   they are flagged with "do not run" - unless ADMIXTOOLS is compiled
   on the system, they will not run.

2. The reviewer reminded me that I should be using `tempdir()` for
   accessing the filesystem. However, in all examples/vignettes/tests
   I've been using `tempdir()` by default unless the user specified
   otherwise.

## Test environments

* My local macOS install: R 4.0.1
* Ubuntu 16.04 GNU/Linux (on Travis CI): R 4.0.0
* Windows (on https://win-builder.r-project.org/): R-devel, R 4.0.2

## R CMD check results

There were no ERRORs and no WARNINGs.

There was only one NOTE ("New submission").

## Downstream dependencies

None.
