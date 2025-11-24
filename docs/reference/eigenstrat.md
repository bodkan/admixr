# EIGENSTRAT data constructor

This function creates an instance of the EIGENSTRAT S3 class, which
encapsulates all paths to data files required for an ADMIXTOOLS
analysis.

## Usage

``` r
eigenstrat(prefix = NULL, ind = NULL, snp = NULL, geno = NULL, exclude = NULL)
```

## Arguments

- prefix:

  Shared path to an EIGENSTRAT trio (set of ind/snp/geno files).

- ind, snp, geno:

  Paths to individual EIGENSTRAT components.

- exclude:

  Pre-defined snp file with excluded sites.

## Value

S3 object of the EIGENSTRAT class.

## Examples

``` r
if (FALSE) # download an example genomic data and get the path prefix to the
# trio of snp/geno/ind files in an EIGENSTRAT format
prefix <- download_data(dirname = tempdir())

# wrap the trio of snp/geno/ind files in an object of the class
# EIGENSTRAT
snps <- eigenstrat(prefix)
#> Error: object 'prefix' not found
 # \dontrun{}
```
