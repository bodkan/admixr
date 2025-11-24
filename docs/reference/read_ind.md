# Read an EIGENSTRAT ind/snp/geno file.

These functions each read one part of the EIGENSTRAT dataset trio.

## Usage

``` r
read_ind(data)

read_snp(data, exclude = FALSE)

read_geno(data)
```

## Arguments

- data:

  EIGENSTRAT data object.

- exclude:

  Read the list of excluded SNPs?

## Value

A data.frame object.

## Details

Note that `read_geno()` will only read plain-text geno files, not
compressed ones.
