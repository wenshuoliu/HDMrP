###-----------missing data in meta-analysis----------###
###Latest edit date: 03/01/2017
###Author: Yajuan Si
setwd("/Users/Shared/ysi//boxsync//Box Sync/projects/meta-analysis/code/")
########################################################
library()

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
#"disttiming", "loctiming","fu_stat"))
#variables subject to missing values

data0$age_cat<-2
data0$age_cat[data$age<50]<-1
data0$age_cat[data$age>=70]<-3


#outcome:
table(alliance$yesdied,useNA="always")
#0: alive 1: dead 2:unknown NAN
table(alliance$LocalRecurrence,useNA="always")
#0: no
table(alliance$DistRecurrence,useNA="always")
#0: no
#failure_time_dist: (distant recurrence time) or (end of follow-up time when no recurrence)
#failure_time_local: (local recurrence time) or (end of follow-up time when no recurrence)
#fu_time: years from study registration to end of follow-up

#covariate

#tomor size
#tumor_size_cat_imp
#1: <2cm 2:2-5cm 3: >5cm or diffuse inflammatory Stage III

# Nodal Status  
#node_pos_imp
# Negative	0
# 1-3 positive nodes	1
# 4-9 positive nodes	2
# >9 positive nodes	3
# Uncertain/Missing	4

# Molecular Subtype Risk Group  
#risk_group_nonodes
# Missing/Unknown	5
# ER or PR +, HER2neu -	1
# ER & PR -, HER2neu-	2
# ER or PR +, HER2neu +	3
# ER & PR -, HER2neu +	4

# Race  race_group
# White	1
# Black or African American	2
# Asian	3
# Native Hawaiian or Other Pacific Islander	4
# American Indian or Alaska Native	5
# Unknown/Other	6

# Ethnicity	ethnicity
# Hispanic or Latino	1
# Not Hispanic or Latino	2
# Unknown	3
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

alpha <- matrix(1,S,K); 
for (s in 1:S){
alpha[s,1:(K-1)]<- seq(-s*0.1,0.5*s,length=K-1)
  #coefficients for studies
}

beta <- matrix(1, p3, K); #coefficients for collected variable
for (j in 1:p3){
  beta[j,1:(K-1)]<-seq(j,-0.5*j,length=K-1) 
}

#latent class generation

c_prob <- exp(Z %*% beta+ S_i %*% alpha)/apply(exp(Z %*% beta + S_i %*% alpha),1,sum)

C_index<-rep(1,N)
for (i in 1:N){
  C_index[i] <- sample(1:K, 1, replace=T, prob=c_prob[i,])
}

#base distribution 

p2<-1; #one missing X
X<-matrix(0,N,p2)
d_j2<-rep(4,p2); #4 levels
psi<-array(0,c(p2,d_j2,K))

for (j2 in 1:p2){
  for (k in 1:K){
    psi[j2,1:(d_j2-1),k]<-2^(-k-1)
    psi[j2,d_j2,k]<- 1-sum(psi[j2,1:(d_j2-1),k])
  }
}

for (i in 1:N){
  for (j2 in 1:p2){
    X[i,j2]<-sample(1:d_j2[j2],1,replace=F,prob=psi[j2,,C_index[i]])
  }
}
X<-as.data.frame(as.factor(X))
names(X)<-"X"
X_matrix<-model.matrix(~X,X)

p1<-2; #two outcome variables  
Y<-matrix(0,N,p1)

#coefficents
theta<-array(0,c(p1,d_j2+p3,K))
for (j1 in 1:p1){
  for (k in 1:K){
    theta[j1,,k]<-0.5*k+j1
  }
}

#variance
sigma<-matrix(1,p1,K)

#Y
for (j1 in 1:p1){
  for (i in 1:N){
    Y[i,j1]<-rnorm(1)* sigma[j1,C_index[i]] + X_matrix[i,]%*%theta[j1,1:d_j2,C_index[i]] + Z[i,]%*%theta[j1,d_j2+1:p3,C_index[i]]
  }
}
#----MCMC----#
#-global parameters-#
nrun <- 10; burn <-0; thin<-1; eff.n<- (nrun-burn)/thin;

#-hyperparameters-#

#-outfiles-#
#for C
beta <- array(0,c(eff.n, p3, K));
alpha <- array(0,c(eff.n, S, K);
#for X
psi<-array(0,c(eff.n,p2,d_j2,K))
#for Y
theta<-array(0,c(eff.n, p1, d_j2+p3, K))
sigma<-array(1,c(eff.n, p1, K))
#-temporary variables-#
pi_k_x <- rep(1,N)
pi_k_y <- rep(1,N)
#-initial values-#

#-sampler-#
for (t in 1:nrun){

  #1) Update latent class indicator
  for (i in 1:N){
    for (k in 1:k){
      pi_k[i]<-exp(Z[i,] %*% beta[t,,k]+ S_i[i,] %*% alpha[t,,k])
  
     pi_k_x[i]<-1
      for (j2 in 1:p2){
        pi_k_x[i]<-pi_k_x[i] * psi[t,j2,X[i,j2],k]
      }
  
      pi_k_y[i]<-1
      for (j1 in 1:p1){
        pi_k_y[i]<-pi_k_y[i] * dnorm(Y[i,j1],X_matrix[i,]%*%theta[t,j1,1:d_j2,k]+ Z[i,]%*%theta[t,j1,d_j2+1:p3,k],sigma[t,j1,k])
      }
  
    pi_k * pi_k_x * pi_k_y
  
  }
  }
}
###------------###
