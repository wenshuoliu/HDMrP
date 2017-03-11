###-----------missing data in meta-analysis----------###
###Latest edit date: 03/09/2017
###Author: Yajuan Si
setwd("/Users/Shared/ysi//boxsync//Box Sync/projects/meta-analysis/code/")
########################################################
library(MASS); library(MCMCpack)

set.seed(20170301)

###------data processing------###
alliance<-readRDS("data/data.rds")
# data0<-subset(alliance, age>0 & race>0, select=c(study,
#                                node_pos_imp, tumor_size_cat_imp,
#                                risk_group_nonodes
#                                #pay_met,er_pr,her2_neu_cat,
#                                age,race,
#                                yesdied,LocalRecurrence,DistRecurrence,
#                                ttrecurrence_loc,ttrecurrence_dist,ttrecurrence_overall,
#                                failure_time_dist,failure_time_local,fu_time)) 
#                               #"disttiming", "loctiming","fu_stat"))
data0<-subset(alliance, select=c(study,node_pos_imp, tumor_size_cat_imp,
                                                 risk_group_nonodes,
                                                 #pay_met,er_pr,her2_neu_cat,
                                                 #age,race,
                                                 yesdied,LocalRecurrence,DistRecurrence,
                                                 ttrecurrence_loc,ttrecurrence_dist,ttrecurrence_overall,
                                                 failure_time_dist,failure_time_local,fu_time)) 
# #"disttiming", "loctiming","fu_stat"))
# #variables subject to missing values
# 
# data0$age_cat<-2
# data0$age_cat[data$age<50]<-1
# data0$age_cat[data$age>=70]<-3
# 
# 
# #outcome:
# table(alliance$yesdied,useNA="always")
# #0: alive 1: dead 2:unknown NAN
# table(alliance$LocalRecurrence,useNA="always")
# #0: no
# table(alliance$DistRecurrence,useNA="always")
# #0: no
# #failure_time_dist: (distant recurrence time) or (end of follow-up time when no recurrence)
# #failure_time_local: (local recurrence time) or (end of follow-up time when no recurrence)
# #fu_time: years from study registration to end of follow-up
# 
# #covariate
# 
# #tomor size
# #tumor_size_cat_imp
# #1: <2cm 2:2-5cm 3: >5cm or diffuse inflammatory Stage III
# 
# # Nodal Status  
# #node_pos_imp
# # Negative	0
# # 1-3 positive nodes	1
# # 4-9 positive nodes	2
# # >9 positive nodes	3
# # Uncertain/Missing	4
# 
# # Molecular Subtype Risk Group  
# #risk_group_nonodes
# # Missing/Unknown	5
# # ER or PR +, HER2neu -	1
# # ER & PR -, HER2neu-	2
# # ER or PR +, HER2neu +	3
# # ER & PR -, HER2neu +	4
# 
# # Race  race_group
# # White	1
# # Black or African American	2
# # Asian	3
# # Native Hawaiian or Other Pacific Islander	4
# # American Indian or Alaska Native	5
# # Unknown/Other	6
# 
# # Ethnicity	ethnicity
# # Hispanic or Latino	1
# # Not Hispanic or Latino	2
# # Unknown	3
###------DPMPM---------###


###-----HDMrP-------###

#----data processing----#

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

alpha_0 <- matrix(1,S,K); 
for (s in 1:S){
alpha_0[s,1:(K-1)]<- seq(-s*0.1,0.5*s,length=K-1)
  #coefficients for studies
}

beta_0 <- matrix(1, p3, K); #coefficients for collected variable
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
#missing data:simulation
miss_ind_Y<-matrix(0,N,p1)
miss_ind_X<-matrix(0,N,p2)
for (j1 in 1:p1){
  miss_ind_Y[,j1]<-(runif(N) <= 0.1)
}
for (j2 in 1:p2){
  miss_ind_X[,j2]<-(runif(N) <= 0.1)
}

#----MCMC----#
#-define function-#
pi_k_fcn<-function(beta,alpha){
  exp(Z %*% beta+ S_i %*% alpha)/apply(exp(Z %*% beta+ S_i %*% alpha),1,sum)
}

#-global parameters-#
nrun <- 1; burn <-0; thin<-1; eff.n<- (nrun-burn)/thin;

#-hyperparameters-#
sigma_theta_0<-1; #prior for theta

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

#use true data without missing
#X<-X
#X_matrix<-model.matrix(~X,X)
#Y<-Y

#-sampler-#
for (t in 1:nrun){

  #1) Update latent class indicator
  for (k in 1:k){
    for (i in 1:N){
      pi_k<-exp(Z %*% beta[t,,k]+ S_i %*% alpha[t,,k])
  
     pi_k_x[i]<-1
      for (j2 in 1:p2){
        pi_k_x[i]<-pi_k_x[i] * psi[t,j2,X[i,j2],k]
      }
  
      pi_k_y[i]<-1
      for (j1 in 1:p1){
        pi_k_y[i]<-pi_k_y[i] * dnorm(Y[i,j1],X_matrix[i,]%*%theta[t,j1,1:d_j2,k]+ Z[i,]%*%theta[t,j1,d_j2+1:p3,k],sigma[t,j1,k])
      }
  
    C_prob[t,,k]<-pi_k * pi_k_x * pi_k_y
    }
    C_prob[t,,] <- C_prob[t,,]/apply(C_prob[t,,],1,sum)
  }

  for (i in 1:N){
    C_i[i]<-sample(1:K,1,replace=F,C_prob[t,i,])
  }
  
 #2) update theta_k
 X_matrix<-model.matrix(~X,X)
 
 X_Z<-as.matrix(cbind(X_matrix,Z))
 
 for (j1 in 1:p1){
   for (k in 1:K){
    X_Z_k<-X_Z[C_i==k,]
    Sigma_theta_k <- solve(t(X_Z_k) %*% X_Z_k + sigma_theta_0 * diag(d_j2+p3))
    Mu_theta_k <- Sigma_theta_k %*% (t(X_Z_k) %*% Y[C_i==k,j1])
    theta[t,j1,,k]<-mvrnorm(1,Mu_theta_k,Sigma_theta_k)
    sigma[t,j1,k]<-sqrt(1/rgamma(1,sum(C_i==k)/2+1,
                            t(Y[C_i==k,j1]-X_Z_k %*% theta[t,j1,,k])%*%(Y[C_i==k,j1]-X_Z_k %*% theta[t,j1,,k])/2+1))
   }
 }
 
 #3) update psi_k, 

 for (k in 1:K){
   for (j2 in 1:p2){
     psi_prob<-rep(0,d_j2[j2]); #temporary
     for (c_x in 1:d_j2[j2]){
       psi_prob[c_x]<-sum(X[C_i==k,j2]==c_x)+1;
     }
     psi_temp<-rgamma(d_j2[j2],psi_prob,1);
     psi[t,j2,,k]<-psi_temp/sum(psi_temp);
     #psi[t,j2,,k]<-rdirichlet(1,psi_prob);
   }
 }
 #4) uodate alpha_k beta_k: noninformative prior
 
 for (k in 1:K){
   if (t==1){
     alpha_new <- mvrnorm(1,alpha_0[,k],diag(S)) #proposal
     beta_new <- mvrnorm(1,beta_0[,k],diag(p3)) 
   }else{
   alpha_new <- mvrnorm(1,alpha[t-1,,k],diag(S)) #proposal
   beta_new <- mvrnorm(1,beta[t-1,,k],diag(p3))
   }
   
   prob_sum_nok<-0
   
   for (h in 1:K){
     if (h !=k){
       if (t==1){
      prob_sum_nok <- prob_sum_nok + pi_k_fcn(beta_0[,h],alpha_0[,h])
       }else{
       prob_sum_nok <- prob_sum_nok + pi_k_fcn(beta[t-1,,h],alpha[t-1,,h])
       }
     }
   }
   
   if (t==1){
    alpha_prob_ratio<-pi_k_fcn(beta_0,alpha_new)^sum(C_i==k)*
       (1-prob_sum_nok-pi_k_fcn(beta_0,alpha_new))^sum(C_i==K)/
       (pi_k_fcn(beta_0,alpha_0)^sum(C_i==k)*)/
       pi_k_fcn(beta_0,alpha_0)^sum(C_i==K)    
  
    
   }else{
   alpha_prob_ratio<-pi_k_fcn(beta[t-1,,k],alpha_new)^sum(C_i==k)*
     (1-prob_sum_nok-pi_k_fcn(beta[t-1,,k],alpha_new))^sum(C_i==K)/
     (pi_k_fcn(beta[t-1,,k],alpha[t-1,,k])^sum(C_i==k)*)/
     pi_k_fcn(beta[t-1,,k],alpha[t-1,,k])^sum(C_i==K)
   
   }
   if (runif(1)< min(alpha_prob_ratio,1)){
     alpha[t,,k]<-alpha_new
   }else{
     alpha[t,,k]<-alpha[t-1,,k]
   }

   if (t==1){
     
   }else{
   beta_prob_ratio<-pi_k_fcn(beta_new, alpha[t,,k])^sum(C_i==k)*
     (1-prob_sum_nok-pi_k_fcn(beta_new, alpha[t,,k]))^sum(C_i==K)/
     (pi_k_fcn(beta[t-1,,k],alpha[t,,k])^sum(C_i==k)*)/
     pi_k_fcn(beta[t-1,,k],alpha[t,,k])^sum(C_i==K)
   }
   
   
   if (runif(1)< min(beta_prob_ratio,1)){
     beta[t,,k]<-beta_new
   }else{
     beta[t,,k]<-beta[t-1,,k]
   }
   
 }

 #5) impute missing Y

 for (j1 in 1:p1){
   miss_index<-(1:N)[miss_ind_Y[,j1]==1] 
   for (i in miss_index){
    Y[i,j1]<-rnorm(1) * sigma[t,j1,C_i[i]] +
     X_matrix[i,]%*%theta[t,j1,1:d_j2,C_i[i]]+ Z[i,]%*%theta[t,j1,d_j2+1:p3,C_i[i]]
   }
 }

 #6) impute missing X
 for (j2 in 1:p2){
   miss_index<-(1:N)[miss_ind_X[,j2]==1]
   for (i in miss_index){
     X[i,j2]<-sample(1:d_j2[j2],1,replace=F,psi[t,j2,,C_i[i]])
   }
 }
}#end of mcmc
###------------###
