# Calculate the D, f4, f4-ratio, or f3 statistic.

Calculate the D, f4, f4-ratio, or f3 statistic.

## Usage

``` r
f4ratio(data, X, A, B, C, O, outdir = NULL, params = NULL)

d(
  data,
  W,
  X,
  Y,
  Z,
  quartets = NULL,
  outdir = NULL,
  f4mode = FALSE,
  params = NULL
)

f4(data, W, X, Y, Z, quartets = NULL, outdir = NULL, params = NULL)

f3(data, A, B, C, outdir = NULL, inbreed = FALSE, params = NULL)
```

## Arguments

- data:

  EIGENSTRAT data object.

- outdir:

  Where to put all generated files (temporary directory by default).

- params:

  Named list of parameters and their values. For instance,
  `params = list(allsnps = "YES")` or `params = list(blgsize = 0.01)`
  (or an arbitrary combination of parameters using a list with multiple
  named elements).

- W, X, Y, Z, A, B, C, O:

  Population names according to the nomenclature used in Patterson et
  al., 2012.

- quartets:

  List of character vectors (quartets of population/sample labels)

- f4mode:

  Calculate the f4 statistic instead of the D statistic.

- inbreed:

  See README.3PopTest in ADMIXTOOLS for an explanation.

## Value

Data frame object with calculated statistics

## Examples

``` r
if (FALSE) # download an example genomic data set and prepare it for analysis
snps <- eigenstrat(download_data(dirname = tempdir()))

# define a set of populations to analyze
pops <- c("French", "Sardinian", "Han", "Papuan", "Dinka")

result_f4ratio <- f4ratio(
    X = pops, A = "Altai", B = "Vindija", C = "Yoruba", O = "Chimp",
    data = snps
)
#> Error: object 'snps' not found

result_d <- d(
    W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp",
    data = snps
)
#> Error: object 'snps' not found

result_f4 <- f4(
    W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp",
    data = snps
)
#> Error: object 'snps' not found

result_f3 <- f3(
    A = pops, B = "Mbuti", C = "Khomani_San",
    data = snps
)
#> Error: object 'snps' not found
 # \dontrun{}
```
