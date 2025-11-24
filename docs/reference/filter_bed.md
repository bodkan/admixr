# Filter EIGENSTRAT data based on a given BED file

Keep (or discard) SNPs that overlap (or lie outside of) regions in a
given BED file.

## Usage

``` r
filter_bed(
  data,
  bed,
  remove = FALSE,
  outfile = tempfile(fileext = ".snp"),
  bedtools_args = ""
)
```

## Arguments

- data:

  EIGENSTRAT data object.

- bed:

  Path to a BED file.

- remove:

  Remove sites falling inside the BED file regions? By default, sites
  that do not overlap BED regions are removed.

- outfile:

  Path to an output snp file with coordinates of excluded sites.

- bedtools_args:

  Optional arguments to \`bedtools intersect\` such as `"-sorted"` or
  `"-sorted -nonamecheck"`.

## Value

Updated S3 EIGENSTRAT data object.

## Details

This function requires a functioning bedtools installation! See:

\- https://github.com/arq5x/bedtools2

\- https://bedtools.readthedocs.io/

## Examples

``` r
if (FALSE) # download an example genomic data set
prefix <- download_data(dirname = tempdir())
# create an EIGENSTRAT R object from the downloaded data
snps <- eigenstrat(prefix)
#> Error: object 'prefix' not found

# get the path to an example BED file
bed <- file.path(dirname(prefix), "regions.bed")
#> Error: object 'prefix' not found

# BED file contains regions to keep in an analysis
snps_kept <- filter_bed(snps, bed)
#> Error: object 'snps' not found
# BED file contains regions to remove from an analysis
snps_removed <- filter_bed(snps, bed, remove = TRUE)
#> Error: object 'snps' not found
 # \dontrun{}
```
