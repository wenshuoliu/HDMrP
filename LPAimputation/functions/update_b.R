###----------------------###
# update random effects in the outcome model #
###----------------------###
update.b <- function(sigma_b=sigma_b, sigma=sigma, T_i = T_i, Y = Y, D_f_rep = D_f_rep, D_r_rep = D_r_rep, t = t, 
                     beta = beta, b_rep = b_rep, alpha_rep = alpha_rep) {
  n <- length(T_i)
  V_b <- rep(0, n)  #variance
  mu_b <- rep(0, n)  #mean
  for (i in 1:n) {
    if (i == 1) {
      j_index <- 1:T_i[1]
    } else {
      j_index <- sum(T_i[1:(i - 1)]) + 1:T_i[i]
    }
    
    V_b[i] <- (sigma_b^(-2) + sigma^(-2) * crossprod(D_r_rep[j_index, ], D_r_rep[j_index, ]))^(-1)
    mu_b[i] <- sigma^(-2) * V_b[i] * crossprod(D_r_rep[j_index, ], Y[j_index] - D_f_rep[j_index, ] %*% beta - t[j_index] * 
                                                 alpha_rep[j_index])
    
  }
  return(rnorm(n) * sqrt(V_b) + mu_b)
}
