# Print the full log output of an admixr wrapper to the console.

Print the full log output of an admixr wrapper to the console.

## Usage

``` r
loginfo(x, target = NA, save = FALSE, prefix = NA, dir = ".", suffix = ".txt")
```

## Arguments

- x:

  Output from one of the admixr wrappers (d, f4, qpAdm, ...)

- target:

  A specific log to examine (relevant for multiple target qpAdm runs)

- save:

  Save the log output to a disk?

- prefix:

  Prefix of the output log file(s) (name of the admixr command by
  default)

- dir:

  In which directory to save the log file(s)?

- suffix:

  Suffix of the output log file(s) (".txt" by default)

## Examples

``` r
if (FALSE) # download an example genomic data set and prepare it for analysis
snps <- eigenstrat(download_data(dirname = tempdir()))

# define a set of populations to analyze and calculate a D statistic
pops <- c("French", "Sardinian", "Han", "Papuan", "Khomani_San", "Mbuti", "Dinka")
result_d <- d(
    W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp",
    data = snps
)
#> Error: object 'snps' not found

# examine the full log output associated with the returned object
loginfo(result_d)
#> Error: object 'result_d' not found
 # \dontrun{}
```
