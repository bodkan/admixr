## Resubmission

1. After reading previous review comments, I'm afraid there is some
   misunderstanding in what this R package does.

   As written in the DESCRIPTION, it's only purpose is to wrap around
   a set of command-line utilities distributed in the ADMIXTOOLS
   software suite (ADMIXTOOLS is a large C codebase).

   Naturally, the examples in the documentation require working
   ADMIXTOOLS installation - all documented functions are convenience
   wrappers and automated pipelines for complex command-line
   operations.

   This is why they are flagged with "do not run" - without compiled
   ADMIXTOOLS command-line programs, they will simply not run.

2. The reviewer reminded me that it's not allowed to modify user home
   filespace in examples/vignettes/tests and that I should be using
   `tempdir()`. However, I did not find a case where this is not the
   case. In all examples/vignettes/tests I'm using `tempdir()` unless
   the user specifies otherwise.

## Test environments

* My local macOS install: R 4.0.1
* Ubuntu 16.04 GNU/Linux (on Travis CI): R 4.0.0
* Windows (on https://win-builder.r-project.org/): R-devel, R 4.0.2

## R CMD check results

There were no ERRORs and no WARNINGs.

There was only one NOTE ("New submission").

## Downstream dependencies

None.
