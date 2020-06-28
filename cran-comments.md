## Resubmission

This is a resubmission of my previous attempt to submit the package.

In this version I have addressed both points raised during review:

* I have added examples to the documentation of each exported function.

* The only calls to `cat()` are now in print.* methods and functions
  whose only purpose is to explicitly write output to the console for
  interactive usage. All other communication with the user is done via
  message(), warning() or stop() functions.

## Test environments

* My local macOS install: R 4.0.1
* Ubuntu 16.04 GNU/Linux (on Travis CI): R 4.0.0
* Debian GNU/Linux (on R-hub): R-devel, R 4.0.0
* Windows (on https://win-builder.r-project.org/): R-devel, R 4.0.0

## R CMD check results

There were no ERRORs and no WARNINGs.

There was only one NOTE ("New submission").

## Downstream dependencies

None.
