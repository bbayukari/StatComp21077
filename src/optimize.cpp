#include <Rcpp.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
#include <Eigen/Eigen>

//' @title example
//' @description A Gibbs sampler using Rcpp, the target density has been given.
//' @return Constant matrix.
//' @export
//[[Rcpp::export]]
Eigen::MatrixXd myfun () {
    Eigen::MatrixXd A(3,2);
    A << 1, 2, 3, 4, 5, 6;
    return A;
}
