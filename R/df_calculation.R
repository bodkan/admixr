f4_df <- function(sites, W, X, Y, Z) {
  a <- sample_freq(sites, W)
  b <- sample_freq(sites, X)
  c <- sample_freq(sites, Y)
  d <- sample_freq(sites, Z)
  
  mean((a - b) * (c - d), na.rm = TRUE)
}



f4ratio_df <- function(sites, X, A, B, C, O) {
  tibble::tibble(
    name = X,
    nea = purrr::map_dbl(X, ~ f4_df(sites, A, O, .x, C) / f4_df(sites, A, O, B, C))
  )
}



d_df <- function(sites, W, X, Y, Z) {
  a <- sample_freq(sites, W)
  b <- sample_freq(sites, X)
  c <- sample_freq(sites, Y)
  d <- sample_freq(sites, Z)
  
  nom <- sum((a - b) * (c - d), na.rm = TRUE)
  denom <- sum((a + b - 2 * a * b) * (c + d - 2 * c * d), na.rm = TRUE)
  
  nom / denom
}



sample_freq <- function(sites, samples) {
  base::rowMeans(sites[, samples]) / 2
}