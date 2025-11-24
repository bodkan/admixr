# Merge two sets of EIGENSTRAT datasets

This function utilizes the 'mergeit' command distributed in ADMIXTOOLS.

## Usage

``` r
merge_eigenstrat(merged, a, b, strandcheck = "NO")
```

## Arguments

- merged:

  Prefix of the path to the merged EIGENSTRAT snp/ind/geno trio.

- a, b:

  Two EIGENSTRAT objects to merge.

- strandcheck:

  Deal with potential strand issues? Mostly for historic reasons. For
  details see the README of ADMIXTOOLS convertf.
