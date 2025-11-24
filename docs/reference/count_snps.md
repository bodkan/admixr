# Count the number/proportion of present/missing sites in each sample

Count the number/proportion of present/missing sites in each sample

## Usage

``` r
count_snps(data, missing = FALSE, prop = FALSE)
```

## Arguments

- data:

  EIGENSTRAT data object.

- missing:

  Count present SNPs or missing SNPs?

- prop:

  Calculate the proportion instead of counts?

## Value

A data.frame object with SNP counts/proportions.

## Examples

``` r
if (FALSE) snps <- eigenstrat(download_data(dirname = tempdir()))

present_count <- count_snps(snps)
#> Error: object 'snps' not found
missing_count <- count_snps(snps, missing = TRUE)
#> Error: object 'snps' not found

present_proportion <- count_snps(snps, prop = TRUE)
#> Error: object 'snps' not found
missing_proportion <- count_snps(snps, missing = TRUE, prop = TRUE)
#> Error: object 'snps' not found
 # \dontrun{}
```
