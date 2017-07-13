###----------------------###
# update sd in the outcome model #
###----------------------###

update.sigma <- function(sigma_a = sigma_a, sigma_b = sigma_b, Y = Y, D_f_rep = D_f_rep, D_r_rep = D_r_rep, 
                         t = t, beta = beta, b_rep = b_rep, alpha_rep = alpha_rep) {
  N <- length(Y)
  inv_sigma2 <- rgamma(1, N/2 + sigma_a, crossprod(Y - D_f_rep %*% beta - D_r_rep * b_rep - t * alpha_rep, 
                                                   Y - D_f_rep %*% beta - D_r_rep * b_rep - t * alpha_rep)/2 + sigma_b)
  return(1/sqrt(inv_sigma2))
}