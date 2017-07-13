###----------------------###
# update the sd pf random effects in the outcome model #
###----------------------###
update.sigma_b <- function(sigmab_a = sigmab_a, sigmab_b = sigmab_a, b = b) {
  n <- length(b)
  invsigmab2 <- rgamma(1, n/2 + sigmab_a, crossprod(b, b)/2 + sigmab_b)
  return(1/sqrt(invsigmab2))
}