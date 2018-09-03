# admixr 0.4.0

* The package now has a tutorial vignette describing the main functionality.
* Simple SNP dataset is now included with the package.
* The API of many utility functions has been simplified and their internals
  re-written.
* `filter_sites` is now implemented using `data.table` and allows overlap with
  an arbitrary BED file.

# admixr 0.3.0

* All wrappers have been given simpler names (`qpDstat()` -> `d()`,
  `qpF4ratio()` -> `f4ratio()`, etc).
* F4 statistic can now be calculated using a separate `f4()` function (`f4mode`
  parameter remains in the `d()` function though, as `f4()` calls `d()`
  internally).
* All tests are performed on Travis CI using installed and compiled ADMIXTOOLS
  software.

# admixr 0.2.0

* The package now includes qpAdm functionality.

* Formal tests for all implemented wrapper functions have been implemented.

# admixr 0.1.0

Implemented:

* D statistic as qpDstat,
* f4 statistic as qpDstat(f4mode = TRUE),
* f4-ratio statistic as qpF4Ratio,
* f3 statistic as qp3Pop.
