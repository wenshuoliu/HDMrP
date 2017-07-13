###-----------Simuate incomplete longitudinal data with latent classes---------###
# Author: YS Latest edit date: 07/10/2017
#setwd("/Users/Shared/ysi/boxsync/Box Sync/projects/diabetes-missingEHRS/code")
setwd("~/Box Sync/projects/diabetes-missingEHRS/code")

rm(list = ls())
set.seed(20170524)
###---------------------------------###
# library
library(mvtnorm)
library(BayesLogit)
###---------------------------------###
# defined functions 
source ("functions/replication.R")
# read data
source ("functions/data_input.R")
# defined functions 
source ("functions/replication.R")
# MCMC functions
source ("functions/update_eta.R")
source ("functions/update_C.R")
source ("functions/update_alpha.R")
source ("functions/update_beta.R")
source ("functions/update_b.R")
source ("functions/update_sigma_b.R")
source ("functions/update_sigma_y.R")
###---------------------------------###
# data
D_r_rep <- rep_fcn(D_r, T_i)
D_f_rep <- rep_fcn(D_f, T_i)
n <- dim(D_c)[1]
# hyperameters
L <- 3
R <- 1000 # total #iterations
burn <- 0 # #burn in
thin <- 1 # thin 

# output
Alpha <- matrix(0, (R - burn)/thin, L)
Beta <- matrix(0, (R - burn)/thin, dim(D_f)[2])
Sigma <- matrix(0, (R - burn)/thin, 1)
Eta <- matrix(0, (R - burn)/thin, L*dim(D_c)[2])

Cs <- matrix(0, (R - burn)/thin, n)
B <- matrix(0, (R - burn)/thin, n)
B_s <- matrix(0, (R - burn)/thin, 1)

# initial values
#not sensitive
sigma <- 1 
beta <- matrix(0,dim(D_f)[2],1) 
eta <- matrix(0, dim(D_c)[2], L) 
# for (c in 2:L) {
#   eta[, c] <- -0.5 * c + 2 * (1:dim(D_c)[2])
# }

#sensitive
C<- rep(0,n)
pi_0 <- exp(D_c %*% eta)/apply(exp(D_c %*% eta), 1, sum)
for (i in 1:n){
  C[i] <- sample(1:L,1,replace=F,prob=pi_0[i,])
}

sigma_b <- 1 #ok
b <- rnorm(n) * sigma_b
b_rep <- rep_fcn(b, T_i)

#alpha <- 1:L -1 #c(1,0,-1)

###---------------------------------###
# MCMC
system.time(
for (it in 1:R) {
  # update alpha,for l=1...L
  C_rep <- rep_fcn(C, T_i)
  alpha <- update.alpha(s0_alpha = 100, sigma = sigma, C_rep = C_rep, t = t, Y = Y, D_f_rep = D_f_rep, 
                        beta = beta, D_r_rep = D_r_rep, b_rep = b_rep)
  # update C
  C <- update.C(D_c = D_c, T_i = T_i, Y = Y, D_r_rep = D_r_rep, D_f_rep = D_f_rep, b_rep = b_rep, eta = eta, 
                alpha = alpha, beta=beta, sigma = sigma)
  alpha_rep <- rep_fcn(alpha[C], T_i)
  
  # update eta
  eta <- update.eta(s=100, eta = eta, D_c = D_c, C = C)
   
  # update beta
  beta <- update.beta(s0_beta = 100, sigma = sigma, D_f_rep = D_f_rep, D_r_rep = D_r_rep, t = t, Y = Y)
  
  # update sigma
  sigma <- update.sigma(sigma_a = 0.001, sigma_b = 0.001, Y = Y, D_f_rep = D_f_rep, D_r_rep = D_r_rep, 
                        t = t, beta = beta, b_rep = b_rep, alpha_rep = alpha_rep)
  
  # update b
  b <- update.b(sigma_b=sigma_b, sigma=sigma, T_i = T_i, Y = Y, D_f_rep = D_f_rep, D_r_rep = D_r_rep, t = t, beta = beta, 
                b_rep = b_rep, alpha_rep = alpha_rep)
  b_rep <- rep_fcn(b, T_i) 
  
  # update sigma_b
  sigma_b <- update.sigma_b(sigmab_a = 0.001, sigmab_b = 0.001, b = b)
  
  # Store Results
  if (it > burn & it%%thin == 0) {
    s <- (it - burn)/thin
    
    Alpha[s, ] <- c(t(alpha))
    Beta[s, ] <- c(t(beta))
    Eta[s, ] <- c(as.vector(eta)) #by column
    
    Cs[s, ] <- C
    B[s, ] <- b
    
    Sigma[s] <- sigma
    B_s[s] <- sigma_b
  }
  
  if (it%%500 == 0) {
    print(it)
  }

}

# End MCMC
) # time

# diagnostics
plot(Sigma)
mean(Sigma)
plot(B_s)
mean(B_s)
#plot(B[,i])
hist(apply(B,2,sd))
plot(Beta[,1])
apply(Beta,2,mean)
plot(Alpha[,2])
apply(Alpha,2,mean)

plot(Eta[,dim(D_c)[2]+1])
apply(Eta,2,mean)

table(Cs[s,])
table(C0)

save.image(file="output/sim.RData")
