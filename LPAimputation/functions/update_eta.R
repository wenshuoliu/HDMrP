###-------------------###
# update component specific coefficients in the multimonial logistic regression #
###-------------------###

#Polya-Gamma
update.eta <- function(s=s, eta = eta, D_c = D_c, C = C){
  n <- nrow(D_c)  # Number of subjects
  q3 <- ncol(D_c)  # Number of regression parms
  if (is.vector(eta) == TRUE) 
    eta <- matrix(eta, nrow = q3, ncol = 1) # Convert eta to matrix if vector
  L <- ncol(eta)  # Number of classes
  r <- matrix(0,n,L-1) #PG parameters for w
  w <- matrix(0,n,L-1) #PG variables
  
  for (l in 2:L){
    r[,l-1] <- D_c %*% eta[,l] - log(apply(exp(D_c %*% eta[,-l]),1,sum))
    w[,l-1] <- rpg(n,1,r[,l-1])
    m_l0 <- as.matrix(as.numeric(C==l)-1/2 + w[,l-1] * log(apply(exp(D_c %*% eta[,-l]),1,sum)),n,1)
    S_l <- solve(t(D_c) %*% diag(w[,l-1]) %*% D_c + 1/s * diag(q3))
    m_l <- S_l %*% t(D_c) %*% m_l0
    eta[,l] <- rmvnorm(1,m_l,S_l)
  }
  
  # for (l in 2:L){
  #   r[,l-1] <- D_c %*% eta[,l] - log(apply(exp(D_c %*% eta[,-l]),1,sum))
  #   for (i in 1:n){
  #   w[i,l-1] <- rpg(1,1,r[i,l-1])
  #   }
  #   m_l0 <- as.matrix(as.numeric(C==l)-1/2 + w[,l-1] * log(apply(exp(D_c %*% eta[,-l]),1,sum)),n,1)
  #   S_l <- solve(t(D_c) %*% diag(w[,l-1]) %*% D_c + 1/s * diag(q3))
  #   m_l <- S_l %*% t(D_c) %*% m_l0
  #   eta[,l] <- rmvnorm(1,m_l,S_l)
  # }
  return(eta)
}
#M-H
# update.eta <- function(eta = eta, D_c = D_c, C = C, s = 3/2, burn = burn, it = it, Eta = Eta) {
#   
#   n <- nrow(D_c)  # Number of subjects
#   q3 <- ncol(D_c)  # Number of regression parms
#   if (is.vector(eta) == TRUE) 
#     eta <- matrix(eta, nrow = q3, ncol = 1)
#   # Convert eta to matrix if vector
#   L <- ncol(eta)  # Number of classes
#   mu0 <- rep(0, q3)  # Prior mean of gamma (for all classes)\t
#   S0 <- s * diag(q3)  # Assume diffuse prior variance
#   lold <- lnew <- rep(0, L)  # Old and new likelihoods
#   lpold <- lpnew <- rep(0, L)  # log priors for eta
#   eta_new <- matrix(0, q3, L)  # Candidate eta
#   
#   covl <- diag(q3)  # Proposal covariance (updated after 2*burn)
#   
#   # Draw Candidate # Note: b/c we tune cov of gnew[k,], need loop
#   for (l in 2:L) {
#     #if (it > 2 * burn + 1) 
#     #  covl <- cov(Eta[(burn + 1):(2 * burn), ((l - 1) * q3 + 1):(l * q3)])
#     
#     eta_new[, l] <- eta[, l] + rmvt(1, sigma = 3 * covl, 3) #2.4/sqrt(q3)
#     # Draw from symmetric MV t-distv with covariance matrix
#     lpold[l] <- dmvnorm(eta[, l], mu0, S0, log = T)
#     # Old Log Prior
#     lpnew[l] <- dmvnorm(eta_new[, l], mu0, S0, log = T)
#     # New Log Prior
#   }
#   
#   # Old and New Class Probabilities (from G-Logit Multinomial)
#   pold <- exp(D_c %*% eta)/apply(exp(D_c %*% eta), 1, sum)
#   pnew <- exp(D_c %*% eta_new)/apply(exp(D_c %*% eta_new), 1, sum)
#   
#   # Old and New Multinomial Likelihoods
#   for (l in 1:L) {
#     lold[l] <- sum(log(pold[C == l, l]))
#     # Class l loglike = sum[I(c_i=l)*log(p_il)]
#     lnew[l] <- sum(log(pnew[C == l, l]))
#   }
#   
#   # MH Acceptance Ratio on Log Scale
#   ratio <- sum(lnew) + sum(lpnew) - (sum(lold) + sum(lpold))
#   if (log(runif(1)) < ratio) 
#     eta <- eta_new
#   
#   return(eta)
#   
# }
