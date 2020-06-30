## Resubmission

This is a third (re-)submission of the package to CRAN.

There was one major point raised in the last review: "your examples in
the documentation are all wrapped in 'do not run'".

I should clarify the following:

The purpose of this package is to wrap around command-line programs in
the ADMIXTOOLS suite for population genetics (a large C codebase, only
available for unix/Linux).

All examples require working ADMIXTOOLS installation. This is why they
are flagged with "do not run" - unless ADMIXTOOLS is compiled on the
system, they will not run.
   
This is easily solved on a service such as Travis CI (all functions in
this package are automatically tested with each new version against
the latest compiled ADMIXTOOLS, including all examples and
vignettes). But that's only because I can compile all dependencies on
the server before running the check.

## Test environments

* My local macOS install: R 4.0.1
* Ubuntu 16.04 GNU/Linux (on Travis CI): R 4.0.0
* Windows (on https://win-builder.r-project.org/): R-devel, R 4.0.2

## R CMD check results

There were no ERRORs and no WARNINGs.

There was only one NOTE ("New submission").

## Downstream dependencies

None.
