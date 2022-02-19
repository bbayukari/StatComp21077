#include <Rcpp.h>
#include <RcppEigen.h>

typedef double (*loss_function_ptr)(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta);


// [[Rcpp::export]]
double loss_API(int n,Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::VectorXd beta, SEXP xpsexp) {
    Rcpp::XPtr<loss_function_ptr> xpfun(xpsexp);
    loss_function_ptr loss_function = *xpfun;
    double loss = 0;
    for(int i = 0; i < n; i++){
        loss = loss_function(X,y,beta);
    }
    return loss;
}

typedef bool (*loss_gradient_ptr)(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXd &g);


// [[Rcpp::export]]
Eigen::VectorXd gradient_API(int n,Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, SEXP xpsexp) {
    Rcpp::XPtr<loss_gradient_ptr> xpfun(xpsexp);
    loss_gradient_ptr loss_gradient = *xpfun;
    Eigen::VectorXd g = Eigen::VectorXd::Zero(2);
    for(int i = 0; i < n; i++){
        loss_gradient(X,y,beta,g);
    }
    return g;
}
