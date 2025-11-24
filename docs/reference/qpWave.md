# Find the most likely number of ancestry waves using the qpWave method.

Given a set of 'left' populations, estimate the lowest number of
necessary admixture sources related to the set of 'right' populations.

## Usage

``` r
qpWave(
  data,
  left,
  right,
  maxrank = NULL,
  details = FALSE,
  outdir = NULL,
  params = NULL
)
```

## Arguments

- data:

  EIGENSTRAT data object.

- left, right:

  Character vectors of populations labels.

- maxrank:

  Maximum rank to test for.

- details:

  Return the A, B matrices used in rank calculations?

- outdir:

  Where to put all generated files (temporary directory by default).

- params:

  Named list of parameters and their values. For instance,
  `params = list(allsnps = "YES")` or `params = list(blgsize = 0.01)`
  (or an arbitrary combination of parameters using a list with multiple
  named elements).

## Value

Table of rank test results.

## Details

It has been shown (Reich, Nature 2012 - Reconstructing Native American
population history) that if the 'left' populations are mixtures of N
different sources related to the set of 'right' populations, the rank of
the matrix of the form \\f_4(left_i, left_j; right_k, right_l)\\ will
have a rank N - 1. This function uses the ADMIXTOOLS command qpWave to
find the lowest possible rank of this matrix that is consistent with the
data.

## Examples

``` r
if (FALSE) # download example data set and prepare it for analysis
snps <- eigenstrat(download_data(dirname = tempdir()))

# run the qpWave wrapper (detailed description in the tutorial vignette)
result <- qpWave(
     left = c("French", "Sardinian", "Han"),
     right = c("Altai", "Yoruba", "Mbuti"),
     data = snps
)
#> Error: object 'snps' not found
 # \dontrun{}
```
