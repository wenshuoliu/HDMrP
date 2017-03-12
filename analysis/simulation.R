####Simulation and testing for the updating functions####

library(MASS)

set.seed(20170301)

###-----HDMrP-------###

#----data processing----#
data0 <- readRDS("analysis/data/data0.rds")

#Z:node_pos_imp, tumor_size_cat_imp
#X:tum_grade,er_pr,her2_neu_cat

#borrow Z from Alliance data
Z<-as.matrix(subset(data0,select=c(node_pos_imp, tumor_size_cat_imp)))
S_i<-model.matrix(~study, data0) #study 369901 as the reference level

#hyperparameter specication
K<-4 # #latent classes
p3 <- dim(Z)[2] # #collected covariates
N <- dim(data0)[1] #total sample size

S<-length(unique(data0$study)) # #studies

alpha_0 <- matrix(0,S,K);
for (s in 1:S){
  alpha_0[s,1:(K-1)]<- seq(-s*0.1,0.5*s,length=K-1)
  #coefficients for studies
}

beta_0 <- matrix(0, p3, K); #coefficients for collected variable
for (j in 1:p3){
  beta_0[j,1:(K-1)]<-seq(j,-0.5*j,length=K-1)
}

#latent class generation

c_prob <- exp(Z %*% beta_0+ S_i %*% alpha_0)/apply(exp(Z %*% beta_0 + S_i %*% alpha_0),1,sum)

C_index<-rep(1,N)
for (i in 1:N){
  C_index[i] <- sample(1:K, 1, replace=T, prob=c_prob[i,])
}

#base distribution

p2<-1; #one missing X
X<-matrix(0,N,p2)
d_j2<-rep(4,p2); #4 levels
psi_0<-array(0,c(p2,d_j2,K))

for (j2 in 1:p2){
  for (k in 1:K){
    psi_0[j2,1:(d_j2-1),k]<-2^(-k-1)
    psi_0[j2,d_j2,k]<- 1-sum(psi_0[j2,1:(d_j2-1),k])
  }
}

for (i in 1:N){
  for (j2 in 1:p2){
    X[i,j2]<-sample(1:d_j2[j2],1,replace=F,prob=psi_0[j2,,C_index[i]])
  }
}
X<-as.data.frame(as.factor(X))
names(X)<-"X"

X_matrix<-model.matrix(~X,X)


p1<-2; #two outcome variables
Y<-matrix(0,N,p1)

#coefficents
theta_0<-array(0,c(p1,d_j2+p3,K))
for (j1 in 1:p1){
  for (k in 1:K){
    theta_0[j1,,k]<-0.5*k+j1
  }
}

#variance
sigma_0<-matrix(1,p1,K)

#Y
for (j1 in 1:p1){
  for (i in 1:N){
    Y[i,j1]<-rnorm(1)* sigma_0[j1,C_index[i]] + X_matrix[i,]%*%theta_0[j1,1:d_j2,C_index[i]] + Z[i,]%*%theta_0[j1,d_j2+1:p3,C_index[i]]
  }
}

#----MCMC----#
#-define function-#
pi_k_fcn<-function(beta,alpha){
  exp(Z %*% beta+ S_i %*% alpha)/apply(exp(Z %*% beta+ S_i %*% alpha),1,sum)
}

#-global parameters-#
nrun <- 1; burn <-0; thin<-1; eff.n<- (nrun-burn)/thin;

#-hyperparameters-#
sigma_theta_0<-1 #prior for theta

#-outfiles-#
#for C
beta <- array(0,c(eff.n, p3, K));
alpha <- array(0,c(eff.n, S, K));
C_prob<- array(0,c(eff.n, N, K));
#for X
psi<-array(0,c(eff.n,p2,d_j2,K));
#for Y
theta<-array(0,c(eff.n, p1, d_j2+p3, K));
sigma<-array(1,c(eff.n, p1, K));

#-temporary variables-#
pi_k_x <- rep(1,N);
pi_k_y <- rep(1,N); C_i <- rep(1,N);

#-initial values-#
beta[1,,]<-beta_0
alpha[1,,]<-alpha_0

psi[1,,,]<-psi_0

theta[1,,,]<-theta_0
sigma[1,,]<-sigma_0

#-sampler-#
for (t in 1:nrun){

  #1) Update latent class indicator
  for (k in 1:K){

    pi_k<-exp(Z %*% beta[t,,k]+ S_i %*% alpha[t,,k])

    #-temporary variables-#
    pi_k_x <- rep(1,N);
    pi_k_y <- rep(1,N); C_i <- rep(1,N);

    for (i in 1:N){

      for (j2 in 1:p2){
        pi_k_x[i]<-pi_k_x[i] * psi[t,j2,X[i,j2],k]
      }

      for (j1 in 1:p1){
        pi_k_y[i]<-pi_k_y[i] * dnorm(Y[i,j1],X_matrix[i,]%*%theta[t,j1,1:d_j2,k]+ Z[i,]%*%theta[t,j1,d_j2+1:p3,k],sigma[t,j1,k])
      }

    }

    C_prob[t,,k]<-pi_k * pi_k_x * pi_k_y
  }
  C_prob[t,,] <- C_prob[t,,]/apply(C_prob[t,,],1,sum)


  for (i in 1:N){
    C_i[i]<-sample(1:K,1,replace=F,C_prob[t,i,])
  }
