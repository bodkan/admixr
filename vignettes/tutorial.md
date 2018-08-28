---
title: "admixr - tutorial"
author: "Martin Petr"
date: "2018-08-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---









This vignette describes how to calculate various population admixture
statistics ($D$, $f_4$, etc.) using the `admixr` package.

**A friendly warning**: many of the statistics implemented here can be quite
sensitive to assumptions about population histories or to errors in the data
(causiing spurious correlations between similarly processed samples, aspecially
in case of ancient DNA). If you want to use these methods, you should really
understand the theory behind them! I highly recommend reading Benjamin Peter's
[wonderful
overview](http://www.genetics.org/content/early/2016/02/03/genetics.115.183913)
of the subject and Nick Patterson's [original ADMIXTOOLS
paper](http://www.genetics.org/content/192/3/1065).



## Introduction

[ADMIXTOOLS](http://www.genetics.org/content/192/3/1065) is a widely used
software package for calculating admixture statistics and testing population
admixture hypotheses. However, although powerful and comprehensive, it is not
exactly known for being user-friendly.

A typical ADMIXTOOLS workflow often involves a combination of `sed`/`awk`/shell
scripting and manual editing to create different configuration files. These are
then passed as command-line arguments to one of ADMIXTOOLS commands and control
how to run a particular analysis. The results are then redirected to another
file, which has to be parsed by the user to extract values of interest, often
using command-line utilities again or (worse) by manual copy-pasting.  The
processed results are then analysed in R, Excel or another program.

This workflow is very cumbersome, especially if one wants to explore many
hypotheses involving different combinations of populations. Most importantly,
however, it makes it difficult to follow good practices of reproducible
science, as it is nearly impossible to construct reproducible automated
"pipelines".

This R package makes it possible to perform all stages of an ADMIXTOOLS
analysis entirely from R. It provides a set of convenient functions that
completely remove the need for "low level" configuration of individual
ADMIXTOOLS programs, allowing users to focus on the analysis itself.










## Installation

To install `admixr` from GitHub you need to install the package `devtools`
first. To do this, you can simply run (in R):


```r
install.packages("devtools")
devtools::install_github("bodkan/admixr")
```

Furthermore, if you want to follow the examples in this vignette, you will need
the [tidyverse](https://www.tidyverse.org) collection of packages for
convenient data analysis, which you can install with:


```r
install.packages("tidyverse")
```

When everything is ready, you can run:


```r
library(admixr)
library(tidyverse)
#> ── Attaching packages ────────────────────────────────── tidyverse 1.2.1 ──
#> ✔ ggplot2 3.0.0     ✔ purrr   0.2.5
#> ✔ tibble  1.4.2     ✔ dplyr   0.7.6
#> ✔ tidyr   0.8.1     ✔ stringr 1.3.1
#> ✔ readr   1.1.1     ✔ forcats 0.3.0
#> ── Conflicts ───────────────────────────────────── tidyverse_conflicts() ──
#> ✖ dplyr::filter() masks stats::filter()
#> ✖ dplyr::lag()    masks stats::lag()
```

**Note that in order to use the `admixr` package, you need a working
installation of ADMIXTOOLS!** You can find installation instructions
[here](https://github.com/DReichLab/AdmixTools/blob/master/README.INSTALL).

**Furthermore, you need to make sure that R can find ADMIXTOOLS binaries on the
`$PATH`.** If this is not the case, running `library(admixr)` will show a
warning message with instructions on how to fix this.






## A note about EIGENSTRAT format

ADMIXTOOLS software uses a peculiar set of genetic file formats, which may seem
strange if you are used to working with [VCF
files](http://samtools.github.io/hts-specs/VCFv4.3.pdf).  However, the basic
idea remains the same - we want to store and access SNP data (REF/ALT alleles)
of a set of individuals at a defined set of genomic positions.

EIGENSTRAT data sets always contain three kinds of files:

- `ind` file - specifies the name, sex and population assignment of each sample
- `snp` file - specifies the positions of SNPs, REF/ALT alleles etc.
- `geno` file - contains SNP data (one row per site, one columne per sample) in
  a dense string-based format:
  - 0: individual is homozygous ALT
  - 1: individual is a heterozygote
  - 2: individual is homozygous REF
  - 9: missing data
  
As you can see, a VCF file is essentially a combination of all three files in a
single file. Luckily for us, all three EIGENSTRAT files usually share a common
path and prefix (at least you should try to make it so whenever you work with
them). This allows us to work with just the prefix instead of worrying about
individual files.

As such, all main `admixr` functions accept a `prefix` argument, which
specifies the path and prefix of all three EIGENSTRAT files (you can still work
with individual files if you need to - using `ind`, `snp` and `geno` arguments
of each `admixr` function - but try to avoid that, as it makes your code mode
verbose).

Here is a prefix of a small testing SNP data set that's distributed with
`admixr`. We will be using this data set in the rest of this vignette.


```r
(eigenstrat <- file.path(system.file(package = "admixr", "extdata"), "snps"))
#> [1] "/Users/martin_petr/projects/admixr/inst/extdata/snps"
```

We can verify that there are three files with this prefix, as they should be:


```r
dir(path = dirname(eigenstrat), full.names = TRUE)
#> [1] "/Users/martin_petr/projects/admixr/inst/extdata/snps.geno"
#> [2] "/Users/martin_petr/projects/admixr/inst/extdata/snps.ind" 
#> [3] "/Users/martin_petr/projects/admixr/inst/extdata/snps.snp"
```

Let's look at their contents.

#### `ind` file

```
#> Chimp        U  Chimp
#> Mbuti        U  Mbuti
#> Yoruba       U  Yoruba
#> Khomani_San  U  Khomani_San
#> Han          U  Han
#> Dinka        U  Dinka
#> Sardinian    U  Sardinian
#> Papuan       U  Papuan
#> French       U  French
#> Vindija      U  Vindija
#> Altai        U  Altai
#> Denisova     U  Denisova
```

The first column (sample name) and the third column (population label) are
generally not the same (sample names often have numerical suffixes, etc.), but
be kept them the same for simplicity. Importantly, when specifying
population/sample arguments in `admixr` functions, the information in the third
column is what is used. For example, if you have individuals such as "French1",
"French2", "French3" in the first column of an `ind` file, all three sharing a
"French" population label in the third column, specifying "French" in an
`admixr` command will combine all three samples in a single population, and
will calculate allele frequency from all of them.

#### `snp` file (first 3 lines)

```
#> 1_832756	1	0.008328	832756	T	G
#> 1_838931	1	0.008389	838931	A	C
#> 1_843249	1	0.008432	843249	A	T
```

#### `geno` file (first 3 lines)

```
#> 902021012000
#> 922221211222
#> 922222122222
```










## Philosophy of `admixr`

The goal of `admixr` is to make ADMIXTOOLS analyses as trivial to perform as
possible, without having to worry about par/pop/left/right configuration files
(as they are known in ADMIXTOOLS' jargon) and other low-level details.

The only interface between you and ADMIXTOOLS is the following set of R
functions:

- `d()`
- `f4()`
- `f4ratio()`
- `f3()`
- `qpAdm()`

Anything that would normally require [dozens of lines of shell
scripts](https://gaworkshop.readthedocs.io/en/latest/contents/06_f3/f3.html)
can be often accomplished by running a single line of R code.









The following sections describe the usage of `admixr` on a set of
example analyses that one might be interested in doing.








## $D$ statistic

Let's say we are interested in the following question: _"Which populations
today show evidence of Neanderthal admixture?_

One way of looking at this is using the following D statistic:
$$D(\textrm{present-day human W}, \textrm{African}, \textrm{Neanderthal}, \textrm{Chimp}).$$

All $D$ statistics are based on comparing the proportions of BABA and ABBA
sites patterns observed in data:

$$D = \frac{\textrm{# BABA sites - # ABBA sites}}{\textrm{# BABA sites + # ABBA sites}}.$$

Significant departure of $D$ from zero indicates an excess of allele sharing
between the first and the third population (positive $D$), or an excess of
allele sharing between the second and the third population (negative $D$). If
we get $D$ that is not significantly different from 0, this suggests that the
first and second populations form a clade, and don't differ in their genetic
affinity to the third population.

Therefore, our $D$ statistic above simply tests whether some modern humans
today admixed with Neanderthals, which would increase their genetic affinity to
this archaic group compared to West Africans (whose ancestors never met
Neanderthals).

Let's save the population names first to make the code below more readable:


```r
pops <- c("French", "Sardinian", "Han", "Papuan", "Khomani_San", "Mbuti", "Dinka")
```

Using the `admixr` package we can then calculate the $D$ statistic above simply
by running:


```r
result <- d(W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp", prefix = eigenstrat)
```

The result is a simple `data.frame`:

```r
head(result)
```


|W           |X      |Y       |Z     |       D|   stderr| Zscore|  BABA|  ABBA|  nsnps|
|:-----------|:------|:-------|:-----|-------:|--------:|------:|-----:|-----:|------:|
|French      |Yoruba |Vindija |Chimp |  0.0313| 0.006933|  4.510| 15802| 14844| 487753|
|Sardinian   |Yoruba |Vindija |Chimp |  0.0287| 0.006792|  4.222| 15729| 14852| 487646|
|Han         |Yoruba |Vindija |Chimp |  0.0278| 0.006609|  4.199| 15780| 14928| 487925|
|Papuan      |Yoruba |Vindija |Chimp |  0.0457| 0.006571|  6.953| 16131| 14721| 487694|
|Khomani_San |Yoruba |Vindija |Chimp |  0.0066| 0.006292|  1.051| 16168| 15955| 487564|
|Mbuti       |Yoruba |Vindija |Chimp | -0.0005| 0.006345| -0.074| 15751| 15766| 487642|

We can see that in addition to the input information, this `data.frame`
contains additional columns:

- `D` - $D$ statistic value
- `stderr` - standard error of the $D$ statistic from the block jackknife
- `Zscore` - $Z$-zscore value (number of standard errors the $D$ is from 0,
  i.e. how strongly do we reject the null hypothesis of no admixture)
- `BABA`/`ABBA` - counts of observed site patterns
- `nsnps` - number of SNPs used for the calculation in this row

Output tables from other `admixr` functions follow a very similar format.

While we could certainly make some inferences straight from this table by
looking at the $Z$-scores, tables in general are not the best representation of
this kind of data, especially as the number of samples increases. This is how
we can use the [`ggplot2`](https://ggplot2.tidyverse.org) package to plot the
results:


```r
ggplot(result, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr))
```

![plot of chunk unnamed-chunk-13](figure/unnamed-chunk-13-1.png)

We can see that all three Africans have $D$ values not significantly different
from 0, meaning that the data is consistent with the null hypothesis of no
Neanderthal ancestry in Africans. On the other hand, the test rejects the null
hypothesis for all non-Africans today, suggesting that Neanderthals admixed
with the ancestors of present-day non-Africans. In fact, this is a similar test
to the one that was used as evidence supporting the Neanderthal admixture
hypothesis in the first place!








## $f_4$ statistic

An alternative way of addressing the previous question is to use the $f_4$
statistic, which is very similar to $D$ statistic and can be calculated as:

$$ f_4 = \frac{\textrm{# BABA sites - # ABBA sites}}{\textrm{# sites}}$$
Again, significant departure of $f_4$ from 0 is informative about gene flow,
in an analogous way to $D$ statistic.

To repeat the previous analysis using $f_4$ statistic, we can run:


```r
result <- f4(W = pops, X = "Yoruba", Y = "Vindija", Z = "Chimp", prefix = eigenstrat)
```


```r
head(result)
```


|W           |X      |Y       |Z     |        f4|   stderr| Zscore|  BABA|  ABBA|  nsnps|
|:-----------|:------|:-------|:-----|---------:|--------:|------:|-----:|-----:|------:|
|French      |Yoruba |Vindija |Chimp |  0.001965| 0.000437|  4.501| 15802| 14844| 487753|
|Sardinian   |Yoruba |Vindija |Chimp |  0.001798| 0.000427|  4.209| 15729| 14852| 487646|
|Han         |Yoruba |Vindija |Chimp |  0.001746| 0.000418|  4.178| 15780| 14928| 487925|
|Papuan      |Yoruba |Vindija |Chimp |  0.002890| 0.000417|  6.924| 16131| 14721| 487694|
|Khomani_San |Yoruba |Vindija |Chimp |  0.000436| 0.000415|  1.051| 16168| 15955| 487564|
|Mbuti       |Yoruba |Vindija |Chimp | -0.000030| 0.000410| -0.074| 15751| 15766| 487642|

We can see by comparing this result to the $D$ statistic result above that we
can make the same conclusions.

You might be wondering why we have both $f_4$ and $D$ if they are so similar.
The truth is that $f_4$ is directly informative about the amount of shared
genetic drift (the "branch length") between pairs of populations, which is, in
many cases, a very useful theoretical property. Other than that, it's often a
matter of personal preference and so `admixr` provides a separate functions for
calculating both.







## $f_4$-ratio statistic

Now we know that non-Africans today carry _some_ Neanderthal ancestry. But what if
we want to know _how much_ Neanderthal ancestry they have? What proportion of their
genomes is of Neanderthal origin?

In general, when we are interested in estimating the *proportion* of ancestry
coming some parental lineage, we can use ratio of $f_4$ statistics.

Using the nomenclature of Patterson et al. 2012, we can perform calculate
$f_4$-ratios using the following code (`X` being a vector of samples for which
we want to estimate Neanderthal ancestry):


```r
result <- f4ratio(X = pops, A = "Altai", B = "Vindija", C = "Yoruba", O = "Chimp", prefix = eigenstrat)
```

The ancestry proportion (a number between 0 and 1) is given in the `alpha`
column:


```r
head(result)
```


|A     |B       |X           |C      |O     |    alpha|   stderr| Zscore|
|:-----|:-------|:-----------|:------|:-----|--------:|--------:|------:|
|Altai |Vindija |French      |Yoruba |Chimp | 0.023774| 0.006173|  3.851|
|Altai |Vindija |Sardinian   |Yoruba |Chimp | 0.024468| 0.006079|  4.025|
|Altai |Vindija |Han         |Yoruba |Chimp | 0.022117| 0.005901|  3.748|
|Altai |Vindija |Papuan      |Yoruba |Chimp | 0.037311| 0.005821|  6.410|
|Altai |Vindija |Khomani_San |Yoruba |Chimp | 0.003909| 0.005923|  0.660|
|Altai |Vindija |Mbuti       |Yoruba |Chimp | 0.000319| 0.005721|  0.056|


```r
ggplot(result, aes(fct_reorder(X, alpha), alpha, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_errorbar(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(y = "Neandertal ancestry proportion", x = "present-day individual")
```

![plot of chunk unnamed-chunk-20](figure/unnamed-chunk-20-1.png)

We can make several observations:

- Again, we don't see any significant Neanderthal ancestry in present-day
  Africans (proportion consistent with 0%), which is what we confirmed using
  $D$ and $f_4$ above.
- Present-day non-Africans carry between 2-3% of Neanderthal ancestry.
- We see a much higher proportion of Neanderthal ancestry in people from Papua
  New Guinea - more than 4%!








## $f_3$ statistic

The $f_3$ statistic, also known as the 3-population statistic, is useful
whenever we want to:

1. Estimate the branch length (shared genetic drift) between a pair of
   populations $A$ and $B$ with respect to a common outgroup $C$. In this case,
   the higher the $f_3$ value, the longer the shared evolutionary time between
   $A$ and $B$.
2. Test whether population $C$ is a mixture of two parental populations $A$ and
   $B$. Negative value of the $f_3$ statistic then serves as statistical
   evidence of this admixture.

As an example, imagine we are interested in relative divergence times between
pairs of present-day human populations and want to know in which approximate
order they split of from each other. To address this problem, we could use
$f_3$ statistic by fixing the $C$ outgroup as Chimp, and calculating pairwise
$f_3$ statistics between all pairs of present-day modern humans.



```r
pops <- c("French", "Sardinian", "Han", "Papuan", "Khomani_San", "Mbuti", "Dinka", "Yoruba")

result <- f3(A = pops, B = pops, C = "Chimp", prefix = eigenstrat)
```


```r
head(result)
```


|A      |B           |C     |           f3|       stderr|  Zscore|  nsnps|
|:------|:-----------|:-----|------------:|------------:|-------:|------:|
|French |Sardinian   |Chimp | 1.738202e+15| 4.373915e+12| 397.402| 225492|
|French |Han         |Chimp | 1.658993e+15| 4.153303e+12| 399.440| 229514|
|French |Papuan      |Chimp | 1.626844e+15| 4.323632e+12| 376.268| 226981|
|French |Khomani_San |Chimp | 1.221558e+15| 3.667817e+12| 333.048| 243884|
|French |Mbuti       |Chimp | 1.259897e+15| 3.474916e+12| 362.569| 253162|
|French |Dinka       |Chimp | 1.410976e+15| 3.752041e+12| 376.055| 257987|


```r
# sort the population labels according to an increasing f3 value relative to French
ordered <- filter(result, A == "French") %>% arrange(f3) %>% .[["B"]] %>% c("French")

# plot heatmap of pairwise f3 values
result %>%
  mutate(A = factor(A, levels = ordered),
         B = factor(B, levels = ordered)) %>%
  ggplot(aes(A, B)) + geom_tile(aes(fill = f3))
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24-1.png)

We can see that when we order the heat map labels based on values of pairwise
$f_3$ statistics, the order of population splits pops up beautifully (i.e.
San separated first, followed by Mbuti, etc.).







## qpAdm method

The last ADMIXTOOLS method implemented in `admixr` is qpAdm. Unfortunately, it
is also one that is the most complex and has not been properly described and
peer-reviewed yet. Nevertheless, it seems to have lot of power to disentangle
complex admixture scenarios, and so we included it in our package as well.

Very briefly, qpAdm can be used to estimate admixture proportions coming from a
series of $N$ source *ancestral* populations, assuming we have reference
populations that form clades with those *source* populations that are closer to
them than to any of the specified *outgroup* populations.

Formally, if a Test population has ancestry coming from $N$ ancestral source
populations, with Reference populations being closer to them than are outgroup
populations $O_i$, we can write:

$$f_4(\textrm{Test}, O_a, O_b, O_c) \approx \sum_{i=1}^N \alpha_i f_4(\textrm{Reference}_i, O_a; O_b, O_c),$$

where $\sum_{i=1}^N \alpha_i = 1$ and $\alpha_i \geq 1$ for all $i = 1, ..., N$.

If this looks like black magic to you, I feel your pain and direct you to the
Supplementary Section 9 of [Haak et al.
2015](http://www.nature.com/articles/nature14317), and an [informal
write-up](https://github.com/DReichLab/AdmixTools/blob/master/pdoc.pdf)
distributed with the ADMIXTOOLS software.

Probably the simplest possible case to show that qpAdm works is by returning
to the question of estimating Neanderthal ancestry proportions. Let's define:

- Europeans as the *target* samples to estimate ancestry proportions for
- Vindija Neanderthal and an African as two *reference* populations (two
  potential sources of ancestries in Europeans today)
- *outgroup* populations - Chimp, Altai Neanderthal and Denisovan (which are
  all further from the true ancestral populations - Vindija and African - than
  the *reference* populations)

Assuming all of that, we can run qpAdm with:


```r
result <- qpAdm(
  target = c("French", "Sardinian", "Mbuti", "Dinka"),
  refs = c("Vindija", "Yoruba"),
  outgroups = c("Chimp", "Denisova", "Altai"),
  prefix = eigenstrat
)
```

And we get the result as a `data.frame` again:

```r
result
```


|target    |  nsnps|    Vindija|    Yoruba| stderr_Vindija| stderr_Yoruba|
|:---------|------:|----------:|---------:|--------------:|-------------:|
|French    | 499434|  0.0215749| 0.9784251|          0.006|         0.006|
|Sardinian | 499314|  0.0246924| 0.9753076|          0.006|         0.006|
|Mbuti     | 499334|  0.0009431| 0.9990569|          0.006|         0.006|
|Dinka     | 499362| -0.0027000| 1.0027000|          0.005|         0.005|


Note that we get admixture proportions standard errors for each potential
source in columns. If we compare this result to the $f_4$-ratio values
calculated above, we see that the qpAdm estimates are very close to what we got
earlier.








## Merging populations

What we've been doing so far was calculating statistics for individual samples.
However, it is often useful to treat multiple samples as a single group or
population. `admixr` provides a function called `merge_labels` that does just
that.

Here is an example: let's say we want to run a similar analysis to the one
described in the $D$ statistic section, but we want to treat Europeans,
Africans and archaics as single combined groups. But the `ind` file that we
have does not contain grouped labels - each sample stands on its own:


```
Chimp        U  Chimp
Mbuti        U  Mbuti
Yoruba       U  Yoruba
Khomani_San  U  Khomani_San
Han          U  Han
Dinka        U  Dinka
Sardinian    U  Sardinian
Papuan       U  Papuan
French       U  French
Vindija      U  Vindija
Altai        U  Altai
Denisova     U  Denisova
```

To merge several individual samples under a combined label we can call
`merge_labels` like this:


```r
# paths to the original ind file and a new modified ind file with merged labels
ind_path <- paste0(eigenstrat, ".ind")
modif_path <- paste0(ind_path, ".merged")

merge_labels(
  ind = ind_path,
  modified_ind = modif_path,
  labels = list( # new population labels
    European = c("French", "Sardinian"),
    African = c("Dinka", "Yoruba", "Mbuti", "Khomani_San"),
    Archaic = c("Vindija", "Altai", "Denisova")
  )
)
```

This is what a modified `ind` file generated by `merge_labels` looks like:

```
Chimp        U  Chimp
Mbuti        U  African
Yoruba       U  African
Khomani_San  U  African
Han          U  Han
Dinka        U  African
Sardinian    U  European
Papuan       U  Papuan
French       U  European
Vindija      U  Archaic
Altai        U  Archaic
Denisova     U  Archaic
```

We can then use "European", "African" and "Archaic" labels in any of the
`admixr` wrapper functions described above, we just have to specify a new `ind`
file in addition to the `prefix` argument. For example:


```r
result <- d(W = "European", X = "African", Y = "Archaic", Z = "Chimp",
            prefix = eigenstrat, ind = modif_path)
```

Here is the result, showing (as we've seen above for individual samples) that
Europeans show genetic affinity to archaic humans compared to Africans today:


```r
head(result)
```


|W        |X       |Y       |Z     |      D|   stderr| Zscore|  BABA|  ABBA|  nsnps|
|:--------|:-------|:-------|:-----|------:|--------:|------:|-----:|-----:|------:|
|European |African |Archaic |Chimp | 0.0225| 0.004404|  5.117| 15487| 14805| 489003|

Note that in the `d()` call we provided both path to shared EIGENSTRAT prefix,
but we also specified the path to the modified `ind` file separately. This overides
the unmodified `ind` file that would be normally used, and ADMIXTOOLS picks up
the combined population labels from the modified `ind` file.
