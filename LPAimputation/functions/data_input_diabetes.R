###----------------------###
# data preprocessiong #
###----------------------###
data <- read.csv("~/Box Sync/projects/diabetes-missingEHRS/code/data/newdata/data_timevar_a1cs_output.csv",header = T)

# time invariant: seudo_patient_id c_chf_elix_bf_coh_ind c_chf_rector_bf_coh_ind c_chf_rector_indc_ihd_ccw2010_ind
# time-varying: cohort_number sd_age_val any_event_ind a1c_result_val

plot(a1c$a1c_result_val)

hist(data$a1c_result_val)

source ("~/Box Sync/projects/diabetes-missingEHRS/code/functions/replication.R")
# parameters
n <- 9000  # #patients
T_i <- sample(5:44, n, replace = T)  # #time points
N <- sum(T_i)  #total observations

q1 <- 2  # #covariates with fixed effected
q2 <- 1  #covariates with random effects
Sigma_b <- 1 * diag(q2)  #covariance matrix for random effects

q3 <- 2  #covariates affecting cluster membership

L <- 3  # #latent classes

###---------------------------------###
# predictors

D_f <- matrix(1, n, q1)  #time-invariant covariates with fixed effects

for (j in 2:q1) {  #first column are 1's
  D_f[, j] <- rbinom(n, 1, 0.2 * j)  #binary covariates
}

D_r <- matrix(1, n, q2)  #time-invariant covariates with random effects (random intercepts)

t <- rep(0, N)  # time-varying covariate
t[1:T_i[1]] <- 1:T_i[1] - 1
for (i in 2:n) {
  t[sum(T_i[1:(i - 1)]) + 1:T_i[i]] <- 1:T_i[i] - 1
}

D_c <- matrix(1, n, q3)
for (j in 1:q3) {
  D_c[, j] <- rnorm(n) #rbinom(n, 1, 0.3 * j)  #binary covariates
}

q1 <- 2  # #covariates with fixed effected
q2 <- 1  #covariates with random effects
sigmab_0 <- 2
Sigma_b <- sigmab_0 * diag(q2)  #covariance matrix for random effects

q3 <- 2  #covariates affecting cluster membership

L <- 3  # #latent classes

###---------------------------------###
# predictors

D_f <- matrix(1, n, q1)  #time-invariant covariates with fixed effects

for (j in 2:q1) {  #first column are 1's
  D_f[, j] <- rbinom(n, 1, 0.2 * j)  #binary covariates
}

D_r <- matrix(1, n, q2)  #time-invariant covariates with random effects (random intercepts)

t <- rep(0, N)  # time-varying covariate
t[1:T_i[1]] <- 1:T_i[1] - 1
for (i in 2:n) {
  t[sum(T_i[1:(i - 1)]) + 1:T_i[i]] <- 1:T_i[i] - 1
}

D_c <- matrix(1, n, q3)
for (j in 1:q3) {
  D_c[, j] <- rnorm(n) #rbinom(n, 1, 0.3 * j)  #binary covariates
}

# coefficients
beta <- matrix(2, q1, 1)  #fixed effects

b <- rmvnorm(n, 0, Sigma_b)  #random effects

eta <- matrix(0, q3, L)
for (c in 2:L) {
  eta[, c] <- -0.5 * c + 2 * (1:q3)
}

#alpha <- rep(0, L)  #class-specific coefficients: only 1 time-varying covariate
alpha <- c(0,-2,2)  # -1, 0, 1

# class indicator
c_prob <- exp(D_c %*% eta)/apply(exp(D_c %*% eta), 1, sum)

# c_prob <- c(0.3, 0.1, 0.6)
C0 <- rep(1, n)
for (i in 1:n) {
  C0[i] <- sample(1:L, 1, replace=F,c_prob[i,])
}
# outcome
D_r_rep <- rep_fcn(D_r, T_i)
b_rep <- rep_fcn(b, T_i)
D_f_rep <- rep_fcn(D_f, T_i)
alpha_rep <- rep_fcn(alpha[C0], T_i)

mu <- D_r_rep * b_rep + D_f_rep %*% beta + t * alpha_rep

sigma_y <- 2
Y <- rnorm(N) * sigma_y + mu

rm(list=ls()[! ls() %in% c("T_i","t", "D_f","D_r","D_c","Y","C0")])
