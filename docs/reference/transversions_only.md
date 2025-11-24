# Remove transversions (C-\>T and G-\>A substitutions)

Remove substitutions that are more likely to be a result of ancient DNA
damage (C-\>T and G-\>A substitutions).

## Usage

``` r
transversions_only(data, outfile = tempfile(fileext = ".snp"))
```

## Arguments

- data:

  EIGENSTRAT data object.

- outfile:

  Path to an output snp file with coordinates of excluded sites.

## Value

Updated S3 EIGENSTRAT data object with an additional 'exclude' slot
specifying the path to the set of SNPs to be removed from a downstream
analysis.

## Examples

``` r
if (FALSE) # download an example genomic data set and prepare it for analysis
snps <- eigenstrat(download_data(dirname = tempdir()))

# perform the calculation only on transversions
snps_tv <- transversions_only(snps)
#> Error: object 'snps' not found
results_d <- d(W = "French", X = "Dinka", Y = "Altai", Z = "Chimp", data = snps_tv)
#> Error: object 'snps_tv' not found
 # \dontrun{}
```
