## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("bodkan/admixr")

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("tidyverse")

## ------------------------------------------------------------------------
library(admixr)
library(tidyverse)

## ------------------------------------------------------------------------
(eigen <- file.path(system.file(package = "admixr", "extdata"), "snps"))

## ------------------------------------------------------------------------
dir(path = dirname(eigen))

## ---- echo = FALSE-------------------------------------------------------
cat(system(paste0("head -n 3 ", eigen, ".ind"), intern = TRUE), sep = "\n")

## ---- echo = FALSE-------------------------------------------------------
cat(system(paste0("head -n 3 ", eigen, ".snp"), intern = TRUE), sep = "\n")

## ---- echo = FALSE-------------------------------------------------------
cat(system(paste0("head -n 3 ", eigen, ".geno"), intern = TRUE), sep = "\n")

## ------------------------------------------------------------------------
# place sample names into variables for more readable code
nonafr <- c("French", "Sardinian", "Han", "Papuan")
afr <- c("Khomani_San", "Mbuti", "Yoruba", "Dinka")

## ------------------------------------------------------------------------
df <- d(W = nonafr, X = "Denisova", Y =  afr, Z = "Chimp", prefix = eigen)

## ----eval = FALSE--------------------------------------------------------
#  head(df)

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(head(df))

## ---- fig.width = 7, fig.height = 4--------------------------------------
ggplot(df, aes(fct_reorder(Y, D), D, color = W)) +
  geom_point() +
  labs(x = "African population Y")

## ------------------------------------------------------------------------
df <- f4(W = nonafr, X = "Denisova", Y =  afr, Z = "Chimp", prefix = eigen)

## ----eval = FALSE--------------------------------------------------------
#  head(df)

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(head(df))

## ---- fig.width = 7, fig.height = 4--------------------------------------
ggplot(df, aes(fct_reorder(Y, f4), f4, color = W)) +
  geom_point() +
  labs(x = "African population Y")

## ------------------------------------------------------------------------
df <- f4ratio(
  X = c("French", "Sardinian", "Han", "Papuan", "Dinka", "Mbuti", "Khomani_San"),
  A = "Altai", B = "Vindija", C = "Yoruba", O = "Chimp",
  prefix = eigen
)

## ---- eval=FALSE---------------------------------------------------------
#  df

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(head(df))

## ---- fig.width = 7, fig.height = 4--------------------------------------
ggplot(df, aes(fct_reorder(X, alpha), alpha, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_errorbar(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(y = "Neandertal ancestry proportion", x = "present-day individual")

