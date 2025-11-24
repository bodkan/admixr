# Package index

## Admixture statistics

- [`f4ratio()`](https://bodkan.net/slendr/reference/f4ratio.md)
  [`d()`](https://bodkan.net/slendr/reference/f4ratio.md)
  [`f4()`](https://bodkan.net/slendr/reference/f4ratio.md)
  [`f3()`](https://bodkan.net/slendr/reference/f4ratio.md) : Calculate
  the D, f4, f4-ratio, or f3 statistic.
- [`qpWave()`](https://bodkan.net/slendr/reference/qpWave.md) : Find the
  most likely number of ancestry waves using the qpWave method.
- [`qpAdm()`](https://bodkan.net/slendr/reference/qpAdm.md) : Calculate
  ancestry proportions in a set of target populations.

## qpAdm model exploration

- [`qpAdm_rotation()`](https://bodkan.net/slendr/reference/qpAdm_rotation.md)
  : Fit qpAdm models based on the rotation strategy described in Harney
  et al. 2020 (bioRxiv)
- [`qpAdm_filter()`](https://bodkan.net/slendr/reference/qpAdm_filter.md)
  : Filter qpAdm rotation results for only 'sensible' models

## Examining output log information

- [`loginfo()`](https://bodkan.net/slendr/reference/loginfo.md) : Print
  the full log output of an admixr wrapper to the console.

## Data processing and filtering

- [`eigenstrat()`](https://bodkan.net/slendr/reference/eigenstrat.md) :
  EIGENSTRAT data constructor
- [`filter_bed()`](https://bodkan.net/slendr/reference/filter_bed.md) :
  Filter EIGENSTRAT data based on a given BED file
- [`transversions_only()`](https://bodkan.net/slendr/reference/transversions_only.md)
  : Remove transversions (C-\>T and G-\>A substitutions)
- [`relabel()`](https://bodkan.net/slendr/reference/relabel.md) : Change
  labels of populations or samples
- [`reset()`](https://bodkan.net/slendr/reference/reset.md) : Reset
  modifications to an EIGENSTRAT object
- [`merge_eigenstrat()`](https://bodkan.net/slendr/reference/merge_eigenstrat.md)
  : Merge two sets of EIGENSTRAT datasets
- [`count_snps()`](https://bodkan.net/slendr/reference/count_snps.md) :
  Count the number/proportion of present/missing sites in each sample

## Reading/writing EIGENSTRAT data

- [`read_ind()`](https://bodkan.net/slendr/reference/read_ind.md)
  [`read_snp()`](https://bodkan.net/slendr/reference/read_ind.md)
  [`read_geno()`](https://bodkan.net/slendr/reference/read_ind.md) :
  Read an EIGENSTRAT ind/snp/geno file.
- [`write_ind()`](https://bodkan.net/slendr/reference/write_ind.md)
  [`write_snp()`](https://bodkan.net/slendr/reference/write_ind.md)
  [`write_geno()`](https://bodkan.net/slendr/reference/write_ind.md) :
  Write an EIGENSTRAT ind/snp/geno file.
- [`read_output()`](https://bodkan.net/slendr/reference/read_output.md)
  : Read an output file from one of the ADMIXTOOLS programs.

## Utility functions

- [`download_data()`](https://bodkan.net/slendr/reference/download_data.md)
  : Download example SNP data.
- [`print(`*`<EIGENSTRAT>`*`)`](https://bodkan.net/slendr/reference/print.EIGENSTRAT.md)
  : EIGENSTRAT print method
- [`print(`*`<admixr_result>`*`)`](https://bodkan.net/slendr/reference/print.admixr_result.md)
  : Print out the admixr result object (dataframe or a list) without
  showing the hidden attributes.
