###----------------------###
# update coefficients for time-invariate variables in the outcome model #
###----------------------###

update.beta <- function(s0_beta = s0_beta, sigma = sigma, D_f_rep = D_f_rep, D_r_rep = D_r_rep, t = t, Y = Y) {
  q1 <- dim(D_f_rep)[2]
  V_beta <- solve(s0_beta^(-1) * diag(q1) + sigma^(-2) * crossprod(D_f_rep, D_f_rep))
  mu_beta <- sigma^(-2) * V_beta %*% crossprod(D_f_rep, Y - D_r_rep * b_rep - t * alpha_rep)
  beta <- t(rmvnorm(1, mu_beta, V_beta))
  return(beta)
}