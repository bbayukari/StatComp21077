library(StatComp21077)


x <- matrix(c(1,0,0,2),nrow=2)
beta <- c(10,20)
y <- c(11,41)
code <- "
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::plugins(cpp14)]]

double loss_function(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta) {
    return (X * beta - y).cwiseAbs2().sum();
}

typedef double (*loss_function_ptr)(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta);

// [[Rcpp::export]]
Rcpp::XPtr<loss_function_ptr> get_loss_function() {
    return(Rcpp::XPtr<loss_function_ptr>(new loss_function_ptr(&loss_function)));
}


bool loss_gradient(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXd &g){
    double fx;
    stan::math::gradient([&](auto beta) {
        return (X * beta - y).cwiseAbs2().sum();
    }, beta, fx, g);
    return true;
}

typedef bool (*loss_gradient_ptr)(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &beta, Eigen::VectorXd &g);

// [[Rcpp::export]]
Rcpp::XPtr<loss_gradient_ptr> get_gradient_function() {
    return(Rcpp::XPtr<loss_gradient_ptr>(new loss_gradient_ptr(&loss_gradient)));
}

"

print(demo(x,y,beta,code,1,1))

# 2
# c(-2,-4)
