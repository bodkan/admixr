## Resubmission

This is a third (re-)submission of the package to CRAN.

There were two new points raised in the last review:

1. After reading the first comment ("your examples in the
   documentation are all wrapped in 'do not run'"), I should clarify
   the following:

   The purpose of this package is to wrap around command-line
   utilities in the unix-based ADMIXTOOLS suite for population
   genetics (a large C codebase). The submitted R package allows the
   user to do all the work in R (data processing/statistics/plotting)
   without ever touching bash/awk/sed/perl.

   All examples in the documentation require working ADMIXTOOLS
   installation - the exported functions are wrappers and pipelines
   for more complex command-line operations happening under the hood.
   
   This is why all examples are flagged with "do not run" - they are
   meaningless unless ADMIXTOOLS is compiled on the system and won't
   even run without it.

2. The reviewer reminded me that I should be using `tempdir()` for
   accessing the filesystem. However, in all examples/vignettes/tests
   I've been using `tempdir()` by default unless the user specified
   otherwise. Perhaps this was just a misunderstanding.

## Test environments

* My local macOS install: R 4.0.1
* Ubuntu 16.04 GNU/Linux (on Travis CI): R 4.0.0
* Windows (on https://win-builder.r-project.org/): R-devel, R 4.0.2

## R CMD check results

There were no ERRORs and no WARNINGs.

There was only one NOTE ("New submission").

## Downstream dependencies

None.
