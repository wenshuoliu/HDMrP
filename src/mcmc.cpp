#include<RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;


//update latent class indicators
//'@export
//[[Rcpp::export]]
void update_Ci(int itt, MatrixXd & C_i, MatrixXd & C_prob, const MatrixXd & S_i, const MatrixXd & Z, const MatrixXd & X,
               const MatrixXd Y, const MatrixXd & alpha, const MatrixXd & beta, const MatrixXd theta, const MatrixXd psi){
  Rcout<<"test"<<std::endl;
}
