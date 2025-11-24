# Fit qpAdm models based on the rotation strategy described in Harney et al. 2020 (bioRxiv)

Fit qpAdm models based on the rotation strategy described in Harney et
al. 2020 (bioRxiv)

## Usage

``` r
qpAdm_rotation(
  data,
  target,
  candidates,
  minimize = TRUE,
  nsources = 2,
  ncores = 1,
  fulloutput = FALSE,
  params = NULL
)
```

## Arguments

- data:

  EIGENSTRAT dataset

- target:

  Target population that is modeled as admixed

- candidates:

  Potential candidates for sources and outgroups

- minimize:

  Test also all possible subsets of outgroups? (default TRUE)

- nsources:

  Number of sources to pull from the candidates

- ncores:

  Number of CPU cores to utilize for model fitting

- fulloutput:

  Report also 'ranks' and 'subsets' analysis from qpAdm in addition to
  the admixture proportions results? (default FALSE)

- params:

  Named list of parameters and their values to be passed to
  [`qpAdm()`](https://bodkan.net/slendr/reference/qpAdm.md).

## Value

qpAdm list with proportions, ranks and subsets elements (as with a
traditional qpAdm run) or just the proportions (determined by the value
of the 'fulloutput' argument)

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
 # \dontrun{}
```
