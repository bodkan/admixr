# Calculate ancestry proportions in a set of target populations.

Calculate ancestry proportions in a set of target populations.

## Usage

``` r
qpAdm(
  data,
  target,
  sources,
  outgroups,
  outdir = NULL,
  params = list(allsnps = "YES", summary = "YES", details = "YES")
)
```

## Arguments

- data:

  EIGENSTRAT data object.

- target:

  Vector of target populations (evaluated one at a time).

- sources:

  Source populations related to true ancestors.

- outgroups:

  Outgroup populations.

- outdir:

  Where to put all generated files (temporary directory by default).

- params:

  Named list of parameters and their values. For instance,
  `params = list(allsnps = "YES")` or `params = list(blgsize = 0.01)`
  (or an arbitrary combination of parameters using a list with multiple
  named elements).

## Value

List of three components: 1. estimated ancestry proportions 2. ranks
statistics 3. analysis of patterns (all possible subsets of ancestry
sources).

## Examples

``` r
if (FALSE) # download example data set and prepare it for analysis
snps <- eigenstrat(download_data(dirname = tempdir()))

# estimate the proportion of Neandertal ancestry in a French
# individual and other associated qpAdm statistics (see detailed
# description in the tutorial vignette)
result <- qpAdm(
    target = "French",
    sources = c("Vindija", "Yoruba"),
    outgroups = c("Chimp", "Denisova", "Altai"),
    data = snps
)
#> Error: object 'snps' not found
 # \dontrun{}
```
