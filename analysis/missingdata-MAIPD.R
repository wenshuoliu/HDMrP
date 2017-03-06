###-----------missing data in meta-analysis----------###
###Latest edit date: 03/01/2017
###Author: Yajuan Si
setwd("/Users/Shared/ysi//boxsync//Box Sync/projects/meta-analysis/code/")
########################################################
library()

set.seed(20170301)

###------data processing------###
alliance<-readRDS("data/alliance.rds")
data0<-with(alliance, select=c("study", "node_pos_imp", "tumor_size_imp","age",
                              "disttiming", "loctiming",
                              "er_pr","her2_neu_cat"))
data0$age_cat<-2
data0$age_cat[data$age<50]<-1
data0$age_cat[data$age>=70]<-3

#    value fu_stat_f 1="Alive" 2="Dead"

###------DPMPM---------###


###-----HDMrP-------###

#hyperparameter specication
K<-4 # #latent classes
p3 <- 2 # #collected covariates
N <- 10000 #total sample size
S<-4 # #studies

alpha <- matrix(1,S,K); 
for (s in 1:S-1){
alpha[s,1:(K-1)]<-seq(-2,2,length=K-1) #coefficients for studies
}
beta <- matrix(1, p3, K); #coefficients for collected variable
for (j in 1:p3){
  beta[j,1:(K-1)]<-seq(-j,j,length=K-1) 
}

#collected covariates and study indicators
Z<-matrix(rnorm(N*p3),N,p3)
S_i<-matrix(sample(1:S,N,replace=T),N,1)

#latent class generation

c_prob <- exp(Z %*% beta + S %*% alpha)./colSums(exp(Z %*% beta + S %*% alpha))
###------------###
