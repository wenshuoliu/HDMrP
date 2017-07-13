# replication function
# rewrite time invariant variables as time varying variables
rep_fcn <- function(X,T_i){
  if (is.vector(X) == TRUE) 
    X <- matrix(X, nrow = length(X), ncol = 1)
  n <- length(T_i)
  q <- dim(X)[2]
  X_rep <- matrix(0,sum(T_i),q)
  X_rep[1:T_i[1],] <- matrix(rep(X[1,],T_i[1]),T_i[1],q,byrow = T)
  for (i in 2:n){
    X_rep[sum(T_i[1:(i - 1)]) + 1:T_i[i], ] <- matrix(rep(X[i,],T_i[i]),T_i[i],q,byrow = T)
  }
  return(X_rep)
}