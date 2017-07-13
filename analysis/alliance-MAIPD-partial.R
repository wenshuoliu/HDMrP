###-----------missing data in meta-analysis----------###
###Latest edit date: 03/09/2017
###Author: Yajuan Si
#setwd("/Users/Shared/ysi//boxsync//Box Sync/projects/meta-analysis/code/")
########################################################
library(MASS); #library(MCMCpack)
library(xtable);
library(NPBayesImpute)
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
data0<-subset(alliance, study== "Z0010"|study=="Z0011"|study== "Z1031"|study== "Z1041"|study== "Z1071",
              select=c(study,node_pos_imp, tumor_size_cat_imp,
                                                 risk_group_nonodes,
                                 ttrecurrence_overall,
                                                 #pay_met,er_pr,her2_neu_cat,
                                                 #age,race,
                                                 yesdied,LocalRecurrence,DistRecurrence,
                                                 ttrecurrence_loc,ttrecurrence_dist))
#,ttrecurrence_overall,
                                                 #failure_time_dist,failure_time_local))
                                 #,fu_time)) 
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

#Alliance data  

dp_data<-subset(data0,select=c(study,tumor_size_cat_imp,node_pos_imp,risk_group_nonodes,LocalRecurrence,DistRecurrence))


dp_data[dp_data=='NaN']<-NA
dp_data$node_pos_imp[dp_data$node_pos_imp==4]<-NA
dp_data$risk_group_nonodes[dp_data$risk_group_nonodes==5] <-NA

for (j in 1:dim(dp_data)[2]){
  dp_data[,j]<-droplevels(as.factor(dp_data[,j]))
}

miss_index<-is.na(dp_data)

table(dp_data$study)

model <- CreateModel(dp_data,NULL,5,10000,0.25,0.25)

#run 1 burnins, 2 mcmc iterations and thin every 2 iterations
model$Run(1,20000,50)
#retrieve parameters from the final iteration
result <- model$snapshot
result$k_star
result$nu
table(result$z)
#convert ImputedX matrix to dataframe, using proper factors/names etc.
ImputedX <- GetDataFrame(result$ImputedX,dp_data)
#View(ImputedX)
#Most exhauststic examples can be found in the demo below
#demo(example_short)
#demo(example)
pdf("dpmpm.pdf")
par(mfrow=c(2, 1))
hist(as.numeric(dp_data[,4]),xlab = "",main="Observed",breaks=0:4,ylim=c(0, 4500))

hist(as.numeric(ImputedX[,4]),xlab = "",main="Imputed",breaks=0:4,ylim=c(0, 4500))
dev.off()

data(NYMockexample)

#############
Z<-as.matrix(subset(data0,select=c(tumor_size_cat_imp,LocalRecurrence,DistRecurrence)))

print(xtable(table(data0$study,data0$tumor_size_cat_imp)))
print(xtable(table(data0$study,data0$node_pos_imp)))

S_i<-model.matrix(~study, data0) #study 369901 as the reference level
X_0<-as.matrix(subset(data0,select=c(node_pos_imp,risk_group_nonodes))
               
Y_0<-as.matrix(subset(data0,select=c(LocalRecurrence,DistRecurrence)))

#Y_0<-as.matrix(subset(data0,select=c(ttrecurrence_loc,ttrecurrence_dist)))

#missing values
X<-X_0
X[X_0==5]<-NA
X[X_0=='NaN']<-NA
Y<-Y_0
# Y[Y_0[,1]<0,1]<-NA
# Y[Y_0[,2]<0,2]<-NA
# Y<-log(Y)
#max(Y[,1]) 72.44
Y[Y_0=='NaN']<-NA
print(xtable(table(data0$study,Y[,2])))



#hyperparameter specication
K<-2 # #latent classes
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

# c_prob <- exp(Z %*% beta_0+ S_i %*% alpha_0)/apply(exp(Z %*% beta_0 + S_i %*% alpha_0),1,sum)
# 
# C_index<-rep(1,N)
# for (i in 1:N){
#   C_index[i] <- sample(1:K, 1, replace=T, prob=c_prob[i,])
# }

#base distribution 

p2<-1; #one missing X
#X<-matrix(0,N,p2)

d_j2<-rep(4,p2); #4 levels

# psi_0<-array(0,c(p2,d_j2,K))
# 
# for (j2 in 1:p2){
#   for (k in 1:K){
#     psi_0[j2,1:(d_j2-1),k]<-2^(-k-1)
#     psi_0[j2,d_j2,k]<- 1-sum(psi_0[j2,1:(d_j2-1),k])
#   }
# }

# for (i in 1:N){
#   for (j2 in 1:p2){
#     X[i,j2]<-sample(1:d_j2[j2],1,replace=F,prob=psi_0[j2,,C_index[i]])
#   }
# }
X<-as.data.frame(as.factor(X))
names(X)<-"X"



p1<-2; #two outcome variables  
#Y<-matrix(0,N,p1)

# #coefficents
# theta_0<-array(0,c(p1,d_j2+p3,K))
# for (j1 in 1:p1){
#   for (k in 1:K){
#     theta_0[j1,,k]<-0.5*k+j1
#   }
# }
# 
# #variance
# sigma_0<-matrix(1,p1,K)

#Y
# for (j1 in 1:p1){
#   for (i in 1:N){
#     Y[i,j1]<-rnorm(1)* sigma_0[j1,C_index[i]] + X_matrix[i,]%*%theta_0[j1,1:d_j2,C_index[i]] + Z[i,]%*%theta_0[j1,d_j2+1:p3,C_index[i]]
#   }
# }
##missing data:simulation
# miss_ind_Y<-matrix(0,N,p1)
# miss_ind_X<-matrix(0,N,p2)
# for (j1 in 1:p1){
#   miss_ind_Y[,j1]<-(runif(N) <= 0.1)
# }
# for (j2 in 1:p2){
#   miss_ind_X[,j2]<-(runif(N) <= 0.1)
# }

#----MCMC----#
#-define function-#
pi_k_i_fcn<-function(Z,S_i,beta,alpha){
  exp(Z %*% beta+ S_i %*% alpha)
}

#-global parameters-#
nrun <- 1000; burn <-0; thin<-1; eff.n<- (nrun-burn)/thin;

#-hyperparameters-#
sigma_theta_0<-1; #prior for theta

#-outfiles-#
#for C
beta <- array(0,c(eff.n, p3, K));
alpha <- array(0,c(eff.n, S, K));

C_prob<- array(0,c(eff.n, N, K));
#for X
psi<-array(1/d_j2,c(eff.n,p2,d_j2,K));
#for Y
theta<-array(0,c(eff.n, p1, d_j2+p3, K));
sigma<-array(1,c(eff.n, p1, K));



#-initial values-#
# beta[1,,]<-beta_0
# alpha[1,,]<-alpha_0
# 
# psi[1,,,]<-psi_0
# 
# theta[1,,,]<-theta_0
# sigma[1,,]<-sigma_0

#initial values of missing
miss_ind_X<-is.na(X)
for (j2 in 1:p2){
X[is.na(X[,j2]),j2]<-sample(1:d_j2[j2],sum(is.na(X[,j2])),replace=T,prob=as.numeric(table(X[!is.na(X[,j2]),j2])))
}

X_matrix<-model.matrix(~X,X)

miss_ind_Y<-is.na(Y)

for (j1 in 1:p1){
  Y[is.na(Y[,j1]),j1]<-rnorm(sum(is.na(Y[,j1])))*sd(Y[!is.na(Y[,j1]),j1]) + mean(Y[!is.na(Y[,j1]),j1]) 
}

print(xtable(table(data0$study,miss_ind_Y[,1])))
#use true data without missing
#X<-X
#
#Y<-Y

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
  
 #2) update theta_k

 
 X_Z<-as.matrix(cbind(X_matrix,Z))
 
 for (j1 in 1:p1){
   for (k in 1:K){
    X_Z_k<-X_Z[C_i==k,
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
 
 for (k in 1:(K-1)){
   
   if (t==1){
     alpha_new <- mvrnorm(1,alpha_0[,k],diag(S)) #proposal
     beta_new <- mvrnorm(1,beta_0[,k],diag(p3)) 
   }else{
    alpha_new <- mvrnorm(1,alpha[t-1,,k],diag(S)) #proposal
    beta_new <- mvrnorm(1,beta[t-1,,k],diag(p3))
   }

     alpha_prob_ratio <-1
     
     for (i in 1:N){
          
         prob_sum_nok_alpha<-0 
        for (h in 1:(K-1)){
           if (h !=k){
             if (t==1){
               prob_sum_nok_alpha <- prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,],beta_0[,h],alpha_0[,h])
             }else{
               prob_sum_nok_alpha <- prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,],beta[t-1,,h],alpha[t-1,,h])
             }
           }
        }
         if (t==1){
           alpha_prob_ratio<-alpha_prob_ratio * pi_k_i_fcn(Z[i,], S_i[i,], beta_0[,k],alpha_new)^(C_i[i]==k) /
             (prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,], beta_0[,k],alpha_new))/
             pi_k_i_fcn(Z[i,], S_i[i,], beta_0[,k],alpha_0[,k])^(C_i[i]==k)*
             (prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,], beta_0[,k],alpha_0[,k]))
         }else{
           alpha_prob_ratio<-alpha_prob_ratio * pi_k_i_fcn(Z[i,], S_i[i,], beta[t-1,,k],alpha_new)^(C_i[i]==k) /
             (prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,], beta[t-1,,k],alpha_new))/
             pi_k_i_fcn(Z[i,], S_i[i,], beta[t-1,,k],alpha[t-1,,k])^(C_i[i]==k)*
             (prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,], beta[t-1,,k],alpha[t-1,,k]))
           
         }
     }
     
   if (runif(1)< min(alpha_prob_ratio,1)){
     alpha[t,,k]<-alpha_new
   }else{
     if (t==1){
      alpha[t,,k]<-alpha_0[,k]
     }else{
      alpha[t,,k]<-alpha[t-1,,k]
     }
   }

  #update beta_k 
     beta_prob_ratio <-1
     
     for (i in 1:N){
       
       prob_sum_nok_beta<-0 
       
       for (h in 1:(K-1)){
         if (h !=k){
           if (t==1){
             prob_sum_nok_beta <- prob_sum_nok_beta + pi_k_i_fcn(Z[i,], S_i[i,],beta_0[,h],alpha[t,,h])
           }else{
             prob_sum_nok_beta <- prob_sum_nok_beta + pi_k_i_fcn(Z[i,], S_i[i,],beta[t-1,,h],alpha[t,,h])
           }
         }
       }
       if (t==1){
         beta_prob_ratio<-beta_prob_ratio * pi_k_i_fcn(Z[i,], S_i[i,], beta_new,alpha[t,,k])^(C_i[i]==k) /
           (prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,], beta_new,alpha[t,,k]))/
           pi_k_i_fcn(Z[i,], S_i[i,], beta_0[,k],alpha[t,,k])^(C_i[i]==k)*
           (prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,], beta_0[,k],alpha[t,,k]))
       }else{
         alpha_prob_ratio<-alpha_prob_ratio * pi_k_i_fcn(Z[i,], S_i[i,], beta_new,alpha[t,,k])^(C_i[i]==k) /
           (prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,], beta_new,alpha[t,,k]))/
           pi_k_i_fcn(Z[i,], S_i[i,], beta[t-1,,k],alpha[t,,k])^(C_i[i]==k)*
           (prob_sum_nok_alpha + pi_k_i_fcn(Z[i,], S_i[i,], beta[t-1,,k],alpha[t,,k]))
         
       }
     }
     
     if (runif(1)< min(beta_prob_ratio,1)){
       beta[t,,k]<-beta_new
     }else{
       if (t==1){
         beta[t,,k]<-beta_0[,k]
       }else{
         beta[t,,k]<-beta[t-1,,k]
       }
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
     
     X_matrix<-model.matrix(~X,X)
     
}#end of mcmc
###------------###
