# Calculate f4 statistic from a data.frame.
f4_df <- function(sites, w, x, y, z) {
  a <- sample_freq(sites, w)
  b <- sample_freq(sites, x)
  c <- sample_freq(sites, y)
  d <- sample_freq(sites, z)
  
  mean((a - b) * (c - d), na.rm = TRUE)
}



# Calculate f4-ratio statistic from a data.frame.
f4ratio_df <- function(sites, x, a, b, c, o) {
  tibble::tibble(
    name = x,
    nea = purrr::map_dbl(x, ~ f4_df(sites, a, o, .x, c) / f4_df(sites, a, o, b, c))
  )
}



# Calculate D statistic from a data.frame.
d_df <- function(sites, w, x, y, z) {
  a <- sample_freq(sites, w)
  b <- sample_freq(sites, x)
  c <- sample_freq(sites, y)
  d <- sample_freq(sites, z)
  
  nom <- sum((a - b) * (c - d), na.rm = TRUE)
  denom <- sum((a + b - 2 * a * b) * (c + d - 2 * c * d), na.rm = TRUE)
  
  nom / denom
}



#' Calculate f3 statistics given four sets of samples.
#'
#' @param sites data.table with simulated SNPs
#' @param a Character vector of sample names.
#' @param b Character vector of sample names.
#' @param c Character vector of sample names.
#'
#' @return f3 statistic
#'
#' @export
f3_df <- function(sites, A, B, C) {
  a <- sample_freq(sites, A)
  b <- sample_freq(sites, B)
  c <- sample_freq(sites, C)
  
  mean((a - c) * (b - c), na.rm = TRUE)
}



#' Calculate allele frequencies of a given set of samples at a given set
#' of sites.
#'
#' @param sites data.table with genotypes.
#' @param samples Character vector of sample names.
#'
#' @import data.table
sample_freq <- function(sites, samples) {
  1 - (base::rowMeans(sites[, samples]) / 2)
}