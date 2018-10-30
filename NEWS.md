# admixr 0.6.1

* It turned out that dragging along Rcpp and Boost dependencies just for the
  VCF -> EIGENSTRAT conversion function causes unnecessary complications in
  the installation process. It's not worth having it in the package if it
  would be used only by a small fraction of potential users.

  This function has been removed and the `vcf2eigenstrat` program is maintained
  in its own repository.

# admixr 0.6.0

* Conversion of VCF to EIGENSTRAT format is now implemented in C++ and should
  be approximately infinitely faster than the old conversion function written
  in pure R.
* Conversion of EIGENSTRAT _into_ VCF has been removed.

# admixr 0.5.0

* Added full implementations of `qpWave()` and `qpAdm()` functions.
* `filter_bed()` now implemented simply by calling `bedtools` in the background.
  This turned out to be way faster and memory efficient than the previous
  data.table-based solution.

# admixr 0.4.1

* Fixed missing `group_labels()` update.
* Removed the huge built-in data set. Implemented `download_data()` function
  that fetches the example data set from the web.

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
