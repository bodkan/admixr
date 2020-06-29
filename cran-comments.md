## Resubmission

In this version I have addressed both points raised during the second
round of reviews:

* I have added examples to the documentation of each exported wrapper
  function.

* The only calls to `cat()` are now in `print.*` methods and in
  interactive functions whose only purpose is to print output to the
  console. All other communication with the user is now done via
  `message()`, `warning()` or `stop()` functions.

## Test environments

* My local macOS install: R 4.0.1
* Ubuntu 16.04 GNU/Linux (on Travis CI): R 4.0.0
* Windows (on https://win-builder.r-project.org/): R-devel, R 4.0.2

## R CMD check results

There were no ERRORs and no WARNINGs.

There was only one NOTE ("New submission").

## Downstream dependencies

None.
