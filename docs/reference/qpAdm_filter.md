# Filter qpAdm rotation results for only 'sensible' models

Filter for p-value larger than a specified cuttof and admixture
proportions between 0 and 1.

## Usage

``` r
qpAdm_filter(x, p = 0.05)
```

## Arguments

- x:

  Output of a qpAdm_rotation() function

- p:

  p-value cutoff (default 0: will only filter for sensible admixture
  proportions)

## Value

qpAdm_rotation object filtered down based on p-value

## Examples

``` r
if (FALSE) # download an example genomic data set and prepare it for analysis
snps <- eigenstrat(download_data(dirname = tempdir()))

# find the set of most likely two-source qpAdm models of
# a French individual - produce only the 'proportions'
# qpAdm summary
models <- qpAdm_rotation(
    data = snps,
    target = "French",
    candidates = c("Dinka", "Mbuti", "Yoruba", "Vindija",
                   "Altai", "Denisova", "Chimp"),
    minimize = TRUE,
    nsources = 2,
    ncores = 2,
    fulloutput = FALSE
)
#> Error: object 'snps' not found

# filter out models which can clearly be rejected
fits <- qpAdm_filter(models, p = 0.05)
#> Error: object 'models' not found
 # \dontrun{}
```
