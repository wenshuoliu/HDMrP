###----------------------###
# update latent class indicator #
###----------------------###
#from multinomial distributions

update.C <- function(D_c = D_c, T_i = T_i, Y = Y, D_r_rep = D_r_rep, D_f_rep = D_f_rep, b_rep = b_rep, 
                     eta = eta, alpha = alpha, beta=beta, sigma = sigma) {
  n <- length(T_i)
  C <- rep(0,n)
  N <- length(Y)
  if (is.vector(eta) == TRUE) 
    eta <- matrix(eta, nrow = q3, ncol = 1)  # Convert eta to matrix if vector
  L <- ncol(eta)  # Number of classes
  ld_l <- matrix(0, N, L)  #likelihood evaluated at class l
  C_prob <- matrix(0, n, L)  #posterior probabilities for classes
  
  pi <- exp(D_c %*% eta)/apply(exp(D_c %*% eta), 1, sum)  #class allocation prob
  
  for (l in 1:L) {
    alpha_rep_l <- rep(alpha[l], N)
    ld_l[, l] <- dnorm(Y, D_r_rep * b_rep + D_f_rep %*% beta + t * alpha_rep_l, sigma)
  }
  
  #dmvnorm(Y[1:T_i[1]], mean=D_r_rep[1:T_i[1]] * b_rep[1:T_i[1]] + D_f_rep[1:T_i[1],] %*% beta + t[1:T_i[1]] * alpha_rep_l[1:T_i[1]], sigma=sigma^2*diag(T_i[1]))
  
  for (i in 1:n) {
    if (i==1){
      C_prob[i, ] <- pi[i, ] * apply(ld_l[1:T_i[1], ], 2, prod) 
    }else{
      C_prob[i, ] <- pi[i, ] * apply(ld_l[sum(T_i[1:(i - 1)]) + 1:T_i[i], ], 2, prod)
    }
    if (sum(C_prob[i, ]==0)==L){
      C[i] <- sample(1:L, 1, replace = F)
    }else{
      C[i] <- sample(1:L, 1, replace = F, prob = C_prob[i, ])
    }
  }
  
  return(C)
}
