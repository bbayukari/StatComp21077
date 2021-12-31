#include <Rcpp.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
#include <Eigen/Eigen>


Eigen::VectorXd expect_y(Eigen::MatrixXd &design, Eigen::VectorXd &coef)
{
Eigen::VectorXd eta = design * coef;
//assert(eta.minCoeff() >= 0); // only use expect_y in where this can be guaranteed.
for (int i = 0; i < eta.size(); i++)
{
    if (eta(i) < 1e-20)
    {
    eta(i) = 1e-20;
    }
}
return eta.cwiseInverse(); // EY is E(Y) = g^-1(Xb), where link func g(u)=1/u in Gamma model.
}
Eigen::VectorXd expect_y(Eigen::MatrixXd &data_matrix, Eigen::VectorXd &beta, double &coef0)
{
Eigen::VectorXd eta = data_matrix * beta + Eigen::VectorXd::Ones(data_matrix.rows()) * coef0;
//assert(eta.minCoeff() >= 0); // only use expect_y in where this can be guaranteed.
for (int i = 0; i < eta.size(); i++)
{
    if (eta(i) < 1e-20)
    {
    eta(i) = 1e-20;
    }
}
return eta.cwiseInverse(); // EY is E(Y) = g^-1(Xb), where link func g(u)=1/u in Gamma model.
}
double loss_function(Eigen::MatrixXd &X, Eigen::VectorXd &y, Eigen::VectorXd &weights, Eigen::VectorXd &beta, double &coef0, double lambda)
{
    int n = X.rows();
    Eigen::VectorXd Xbeta = X * beta + Eigen::VectorXd::Ones(n) * coef0;
    for (int i = 0; i < Xbeta.size(); i++)
    {
        if (Xbeta(i) < 1e-20)
        {
            Xbeta(i) = 1e-20;
        }
    }
    return (Xbeta.cwiseProduct(y) - Xbeta.array().log().matrix()).dot(weights) + lambda * beta.cwiseAbs2().sum();
}

//' @title gamma fit approximate newton method
//' @description fit gamma model with approximate newton method.
//' @param x Input matrix, of dimension \eqn{n \times p}; each row is an observation vector and each column is a predictor/feature/variable.
//' @param y The response variable, of \code{n} observations.
//' @param weights Observation weights.
//' @param beta Initial value of linear coefficient.
//' @param coef0 Initial value of intercept.
//' @param lambda L2 penalty coefficient.
//' @return Intercept and linear coefficient.
//[[Rcpp::export]]
Eigen::VectorXd gamma_fit_approximate_newton_method(Eigen::MatrixXd x, Eigen::VectorXd y, Eigen::VectorXd weights, Eigen::VectorXd beta, double coef0, double lambda)
{
    int n = x.rows();
    int p = x.cols();
    int max_iter = 80;
    Eigen::VectorXd coef = Eigen::VectorXd::Zero(p + 1);
    if (p == 0)
    {
        coef(0) = weights.sum() / weights.dot(y);
        return coef;
    }
    // expand data matrix to design matrix
    Eigen::MatrixXd X(n, p + 1);
    X.rightCols(p) = x;
    X.col(0) = Eigen::MatrixXd::Ones(n, 1);
    coef(0) = coef0;
    coef.tail(p) = beta;

    // modify start point to make sure Xbeta > 0
    double min_xb = (X * coef).minCoeff();
    if (min_xb < 1e-20)
    {
        coef(0) += abs(min_xb) + 0.1;
    }

    // Approximate Newton method
    double step = 1;
    Eigen::VectorXd g(p + 1);
    Eigen::VectorXd coef_new;
    Eigen::VectorXd h_diag(p + 1);
    Eigen::VectorXd desend_direction; // coef_new = coef + step * desend_direction
    Eigen::VectorXd EY = expect_y(X, coef);
    Eigen::VectorXd W = EY.array().square() * weights.array();
    coef0 = coef(0);
    beta = coef.tail(p);
    double loglik_new, loglik = -loss_function(x,y,weights,beta,coef0,lambda);

    for (int j = 0; j < max_iter; j++)
    {
        for (int i = 0; i < p + 1; i++)
        {
            h_diag(i) = X.col(i).cwiseProduct(W).dot(X.col(i)) + 2 * lambda; // diag of Hessian
            // we can find h_diag(i) >= 0
            if (h_diag(i) < 1e-7)
            {
                h_diag(i) = 1e7;
            }
            else
                h_diag(i) = 1.0 / h_diag(i);
        }

        g = X.transpose() * (EY - y).cwiseProduct(weights) - 2 * lambda * coef; // negtive gradient direction
        desend_direction = g.cwiseProduct(h_diag);
        coef_new = coef + step * desend_direction; // Approximate Newton method
        coef0 = coef_new(0);
        beta = coef_new.tail(p);
        loglik_new = -loss_function(x,y,weights,beta,coef0,lambda);

        while (loglik_new < loglik && step > 1e-8)
        {
            step = step / 2;
            coef_new = coef + step * desend_direction;
            coef0 = coef_new(0);
            beta = coef_new.tail(p);
            loglik_new = -loss_function(x,y,weights,beta,coef0,lambda);
        }

        bool condition1 = step < 1e-8;
        bool condition2 = abs(loglik - loglik_new) / (0.1 + abs(loglik_new)) < 1e-8;
        bool condition3 = abs(loglik - loglik_new) < 1e-3;
        if (condition1 || condition2 || condition3)
        {
            break;
        }

        coef = coef_new;
        loglik = loglik_new;
        EY = expect_y(X, coef);
        W = EY.array().square() * weights.array();
    }

    return coef;
}


//' @title gamma fit IWLS method
//' @description fit gamma model with IWLS method.
//' @param x Input matrix, of dimension \eqn{n \times p}; each row is an observation vector and each column is a predictor/feature/variable.
//' @param y The response variable, of \code{n} observations.
//' @param weights Observation weights.
//' @param beta Initial value of linear coefficient.
//' @param coef0 Initial value of intercept.
//' @param lambda L2 penalty coefficient.
//' @return Intercept and linear coefficient.
//[[Rcpp::export]]
Eigen::VectorXd gamma_fit_IWLS_method(Eigen::MatrixXd x, Eigen::VectorXd y, Eigen::VectorXd weights, Eigen::VectorXd beta, double coef0, double lambda)
{
    int n = x.rows();
    int p = x.cols();
    int max_iter = 80;
    Eigen::VectorXd coef = Eigen::VectorXd::Zero(p + 1);
    if (p == 0)
    {
        coef(0) = weights.sum() / weights.dot(y);
        return coef;
    }
    // expand data matrix to design matrix
    Eigen::MatrixXd X(n, p + 1);
    X.rightCols(p) = x;
    X.col(0) = Eigen::MatrixXd::Ones(n, 1);
    coef(0) = coef0;
    coef.tail(p) = beta;

    // modify start point to make sure Xbeta > 0
    double min_xb = (X * coef).minCoeff();
    if (min_xb < 1e-20)
    {
        coef(0) += abs(min_xb) + 0.1;
    }


    Eigen::MatrixXd X_new(X);
    Eigen::MatrixXd lambdamat = Eigen::MatrixXd::Identity(p + 1, p + 1);
    lambdamat(0, 0) = 0;
    Eigen::VectorXd EY = expect_y(X, coef);
    Eigen::VectorXd EY_square = EY.array().square();
    Eigen::VectorXd W = EY_square.cwiseProduct(weights); // the weight matriX of IW(eight)LS method
    Eigen::VectorXd Z = X * coef - (y - EY).cwiseQuotient(EY_square);
    coef0 = coef(0);
    beta = coef.tail(p);
    double loglik_new, loglik = -loss_function(x,y,weights,beta,coef0,lambda);

    for (int j = 0; j < max_iter; j++)
    {
        for (int i = 0; i < p + 1; i++)
        {
            X_new.col(i) = X.col(i).cwiseProduct(W);
        }
        Eigen::MatrixXd XTX = 2 * lambda * lambdamat + X_new.transpose() * X;
        coef = XTX.ldlt().solve(X_new.transpose() * Z);
        coef0 = coef(0);
        beta = coef.tail(p);
        loglik_new = -loss_function(x,y,weights,beta,coef0,lambda);
        bool condition2 = abs(loglik - loglik_new) / (0.1 + abs(loglik_new)) < 1e-8;
        bool condition3 = abs(loglik - loglik_new) < 1e-3;
        if (condition2 || condition3)
        {
            break;
        }

        double min_xb = (X * coef).minCoeff();
        if (min_xb < 1e-20)
        {
            coef(0) += abs(min_xb) + 0.1;
        }
        EY = expect_y(X, coef);
        EY_square = EY.array().square();
        loglik = loglik_new;
        W = EY_square.cwiseProduct(weights);
        Z = X * coef - (y - EY).cwiseQuotient(EY_square);
    }

    return coef;
}



