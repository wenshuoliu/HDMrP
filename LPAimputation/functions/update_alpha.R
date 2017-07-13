###----------------------###
# update class specific coefficients in the outcome model #
###----------------------###
update.alpha <- function(s0_alpha = s0_alpha, sigma = sigma, C_rep = C_rep, t = t, Y = Y, D_f_rep = D_f_rep, 
                         beta = beta, D_r_rep = D_r_rep, b_rep = b_rep) {
  
  L <- length(unique(C_rep))
  alpha <- rep(0, L)
  V_alpha <- rep(0, L-1)  #posterior variance for alpha
  mu_alpha <- rep(0, L-1)  #posterior mean for alpha
  for (l in 2:L) {
    V_alpha[l] <- 1/(s0_alpha^(-1) + sigma^(-2) * sum(t[C_rep == l]^2))
    mu_alpha[l] <- V_alpha[l] * sigma^(-2) * sum(t[C_rep == l] * (Y[C_rep == l] - D_f_rep[C_rep == l, ] %*% beta - 
                                                                 D_r_rep[C_rep == l] * b_rep[C_rep == l]))
    alpha[l] <- rnorm(1) * sqrt(V_alpha[l]) + mu_alpha[l] 
  }
  return(alpha)
}
