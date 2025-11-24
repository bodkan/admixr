# Change labels of populations or samples

Replace population/sample names with specified group labels.

## Usage

``` r
relabel(data, ..., outfile = tempfile(fileext = ".ind"))
```

## Arguments

- data:

  EIGENSTRAT trio.

- ...:

  Population/sample names to merge (each new group defined as a
  character vector).

- outfile:

  Path to an output snp file with coordinates of excluded sites.

## Value

Updated S3 EIGENSTRAT data object with an additional 'group' slot
specifying the path to a new ind file. \#'

## Examples

``` r
if (FALSE) # download an example genomic data set and prepare it for analysis
snps <- eigenstrat(download_data(dirname = tempdir()))

# group individual samples into larger populations, creating a new
# EIGENSTRAT R object
new_snps <- relabel(
    snps,
    European = c("French", "Sardinian"),
    African = c("Dinka", "Yoruba", "Mbuti", "Khomani_San"),
    Archaic = c("Vindija", "Altai", "Denisova")
)
#> Error: object 'snps' not found
 # \dontrun{}
```
