#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <stdio.h>
//[[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;
#include <Eigen/Eigen>

//[[Rcpp::export]]
Eigen::MatrixXd pro(Eigen::MatrixXd X,Eigen::MatrixXd beta, Eigen::VectorXd coef0){
  int n = X.rows();
  int p = X.cols();
  int k = coef0.size() - 1;

  Eigen::VectorXd xbeta;
  Eigen::VectorXd coef = Eigen::VectorXd::Zero(p + k);
  coef.head(k) = coef0.head(k);
  coef.tail(p) = beta.col(0);
  Eigen::MatrixXd logit(n, k);
  Eigen::MatrixXd P(n, k + 1);

  xbeta = X * coef.tail(p);
  // compute logit
  for (int i1 = 0; i1 < n; i1++)
  {
    for (int i2 = 0; i2 < k; i2++)
    {
      logit(i1, i2) = 1.0 / (1 + exp(-xbeta(i1) - coef(i2)));
    }
  }
  // compute P
  for (int i1 = 0; i1 < n; i1++)
  {
    for (int i2 = 0; i2 < k + 1; i2++)
    {
      if (i2 == 0)
      {
        P(i1, 0) = logit(i1, 0);
      }
      else if (i2 == k)
      {
        P(i1, k) = 1 - logit(i1, k - 1);
      }
      else
      {
        P(i1, i2) = logit(i1, i2) - logit(i1, i2 - 1);
      }
      if (P(i1, i2) < 1e-10)
        P(i1, i2) = 1e-10;
    }
  }
  return P;
}

//[[Rcpp::export]]
Eigen::MatrixXd grad_cpp(Eigen::MatrixXd X,Eigen::MatrixXd y,Eigen::MatrixXd beta, Eigen::VectorXd coef0){
  int n = X.rows();
  int p = X.cols();
  int k = coef0.size() - 1;

  Eigen::VectorXd xbeta;
  Eigen::VectorXd coef = Eigen::VectorXd::Zero(p + k);
  coef.head(k) = coef0.head(k);
  coef.tail(p) = beta.col(0);
  Eigen::MatrixXd logit(n, k);
  Eigen::MatrixXd P(n, k + 1);
  Eigen::VectorXd g(p + k);
  Eigen::MatrixXd grad_L(n, k);

  xbeta = X * coef.tail(p);
  // compute logit
  for (int i1 = 0; i1 < n; i1++)
  {
    for (int i2 = 0; i2 < k; i2++)
    {
      logit(i1, i2) = 1.0 / (1 + exp(-xbeta(i1) - coef(i2)));
    }
  }
  // compute P
  for (int i1 = 0; i1 < n; i1++)
  {
    for (int i2 = 0; i2 < k + 1; i2++)
    {
      if (i2 == 0)
      {
        P(i1, 0) = logit(i1, 0);
      }
      else if (i2 == k)
      {
        P(i1, k) = 1 - logit(i1, k - 1);
      }
      else
      {
        P(i1, i2) = logit(i1, i2) - logit(i1, i2 - 1);
      }
      if (P(i1, i2) < 1e-10)
        P(i1, i2) = 1e-10;
    }
  }
  // compute gradient

  for (int i1 = 0; i1 < n; i1++)
  {
    for (int i2 = 0; i2 < k; i2++)
    {
      grad_L(i1, i2) = (y(i1, i2) / P(i1, i2) - y(i1, i2 + 1) / P(i1, i2 + 1)) * logit(i1, i2) * (1.0 - logit(i1, i2));
    }
  }
  g.head(k) = grad_L.colwise().sum();
  g.tail(p) = grad_L.rowwise().sum().transpose() * X;

  return g;

}


//[[Rcpp::export]]
double loss_function(Eigen::MatrixXd &X, Eigen::MatrixXd &y, Eigen::VectorXd &weights, Eigen::MatrixXd &beta, Eigen::VectorXd &coef0, double lambda)
{
    int n = X.rows();
    int p = X.cols();
    int k = coef0.size() - 1;

    Eigen::VectorXd xbeta = X * beta.col(0);
    double loss = lambda * beta.col(0).cwiseAbs2().sum();

    double pro = 0;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < k + 1; j++)
      {
        if (y(i, j) == 1)
        {
          if (j == 0)
          {
            loss += log(1 + exp(-xbeta(i) - coef0(0)));
          }
          else if (j == k)
          {
            loss -= log(1 - 1.0 / (1 + exp(-xbeta(i) - coef0(k - 1))));
          }
          else
          {
            pro = 1.0 / (1 + exp(-xbeta(i) - coef0(j))) - 1.0 / (1 + exp(-xbeta(i) - coef0(j - 1)));
            if (pro < 1e-20)
              pro = 1e-20;
            loss -= log(pro);
          }
          break;
        }
      }
    }

    for (int i = 1; i < k; i++)
    {
      if (coef0(i) <= coef0(i - 1))
      {
        cout << "warning: intercept isn't parallel!"
             << " loss is " << loss << endl; //test
        cout << "coef0 are " << coef0 << endl;
        break;
      }
    }

    return loss;
}


//[[Rcpp::export]]
Eigen::VectorXd fit_ordinal(Eigen::MatrixXd X, Eigen::MatrixXd y, Eigen::VectorXd weights, Eigen::MatrixXd beta, Eigen::VectorXd coef0, double lambda, int primary_model_fit_max_iter,double step0)
  {
    //static int num_ = 0; //test
    //int id = ++num_;

    int i;
    int n = X.rows();
    int p = X.cols();
    int k = coef0.size() - 1;

    printf("%d;%d;%d;%d;%d;%d;%d;%d\n", X.rows(), X.cols(), y.rows(), y.cols(), weights.size(), beta.rows(), beta.cols(), coef0.size());
    //cout << X.rows() << ";" << X.cols() << ";" << y.rows() << ";" << y.cols() << ";" << weights.size() << ";" << beta.rows() << ";" << beta.cols() << ";" << coef0.size() << endl;

    //make sure that coef0 is increasing
    for (int i = 1; i < k; i++)
    {
      if (coef0(i) <= coef0(i - 1))
      {
        coef0(i) = coef0(i - 1) + 1;
      }
    }

    double step = step0;
    double loglik_new, loglik;
    Eigen::VectorXd g(p + k);
    Eigen::VectorXd coef_new;
    Eigen::VectorXd h_diag = Eigen::VectorXd::Zero(k + p);
    Eigen::VectorXd xbeta;
    Eigen::MatrixXd logit(n, k);
    Eigen::MatrixXd P(n, k + 1);
    Eigen::VectorXd W(n);
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(n, k);
    Eigen::MatrixXd dL = Eigen::MatrixXd::Zero(n, k);
    Eigen::VectorXd desend_direction; // coef_new = coef + step * desend_direction
    Eigen::MatrixXd grad_L(n, k);
    Eigen::VectorXd coef = Eigen::VectorXd::Zero(p + k);
    coef.head(k) = coef0.head(k);
    coef.tail(p) = beta.col(0);
    loglik = -loss_function(X, y, weights, beta, coef0, lambda);

    printf("0: loss = %f\n", -loglik);

    for (int j = 0; j < primary_model_fit_max_iter; j++)
    {
      xbeta = X * coef.tail(p);
      // compute logit
      for (int i1 = 0; i1 < n; i1++)
      {
        for (int i2 = 0; i2 < k; i2++)
        {
          logit(i1, i2) = 1.0 / (1 + exp(-xbeta(i1) - coef(i2)));
        }
      }
      // compute P
      for (int i1 = 0; i1 < n; i1++)
      {
        for (int i2 = 0; i2 < k + 1; i2++)
        {
          if (i2 == 0)
          {
            P(i1, 0) = logit(i1, 0);
          }
          else if (i2 == k)
          {
            P(i1, k) = 1 - logit(i1, k - 1);
          }
          else
          {
            P(i1, i2) = logit(i1, i2) - logit(i1, i2 - 1);
          }
          if (P(i1, i2) < 1e-10)
            P(i1, i2) = 1e-10;
        }
      }
      // compute gradient

      for (int i1 = 0; i1 < n; i1++)
      {
        for (int i2 = 0; i2 < k; i2++)
        {
          grad_L(i1, i2) = (y(i1, i2) / P(i1, i2) - y(i1, i2 + 1) / P(i1, i2 + 1)) * logit(i1, i2) * (1.0 - logit(i1, i2));
        }
      }
      g.head(k) = grad_L.colwise().sum();
      g.tail(p) = grad_L.rowwise().sum().transpose() * X - 2 * lambda * coef.tail(p).eval();

      for (int i2 = 0; i2 < k; i2++)
      {
        for (int i1 = 0; i1 < n; i1++)
        {
          h_diag(i2) += (1.0 / P(i1, i2) + 1.0 / P(i1, i2 + 1)) * logit(i1, i2) * logit(i1, i2) * (1.0 - logit(i1, i2)) * (1.0 - logit(i1, i2));
        }
        if (h_diag(i2) < 1e-7 && h_diag(i2) >= 0)
          h_diag(i2) = 1e7;
        else if (h_diag(i2) > -1e-7 && h_diag(i2) < 0)
          h_diag(i2) = -1e7;
        else
          h_diag(i2) = 1.0 / h_diag(i2);
      }
      W = Eigen::VectorXd::Zero(n);
      for (int i = 0; i < n; i++)
      {
        for (int i1 = 0; i1 < k; i1++)
        {
          W(i) += (1.0 / P(i, i1) + 1.0 / P(i, i1 + 1)) * logit(i, i1) * logit(i, i1) * (1.0 - logit(i, i1)) * (1.0 - logit(i, i1));
        }
        for (int i1 = 0; i1 < k - 1; i1++)
        {
          W(i) -= 1.0 / P(i, i1 + 1) * logit(i, i1) * logit(i, i1 + 1) * (1.0 - logit(i, i1)) * (1.0 - logit(i, i1 + 1));
        }
      }

      for (int i = 0; i < p; i++)
      {
        h_diag(i + k) = X.col(i).cwiseProduct(W).dot(X.col(i)) + 2 * lambda;
        if (h_diag(i + k) < 1e-7 && h_diag(i + k) >= 0)
          h_diag(i + k) = 1e7;
        else if (h_diag(i + k) > -1e-7 && h_diag(i + k) < 0)
          h_diag(i + k) = -1e7;
        else
          h_diag(i + k) = 1.0 / h_diag(i + k);
      }
      /*
            if(id==5&&j==0){
              cout << "xbeta:" << endl << xbeta << endl;
              cout << "logit:" << endl << logit << endl;
              cout << "P:" << endl << P << endl;
              cout << "grad_L:" << endl << grad_L << endl;
            }
      */
      //step = 1;
      //desend_direction = g;
      desend_direction = g.cwiseProduct(h_diag);
      coef_new = coef + step * desend_direction; // ApproXimate Newton method
      printf("%d: step = %f, coef0 = ",j, step);
      for (int i1 = 0; i1 < k; i1++)
      {
        printf("%f,", coef_new(i1));
      }
      printf("\n");
      i = 1;
      while(i < k)
      {
        for (i = 1; i < k; i++)
        {
          if (coef_new(i) <= coef_new(i - 1))
          {
            step = step / 2;
            coef_new = coef + step * desend_direction;
            break;
          }
        }
      }

      beta.col(0) = coef_new.tail(p);
      coef0.head(k) = coef_new.head(k);
      loglik_new = -loss_function(X, y, weights, beta, coef0, lambda);
      printf("%d: step = %f, coef0 = ", j, step);
      for (int i1 = 0; i1 < k; i1++)
      {
        printf("%f,", coef_new(i1));
      }
      printf("\n");
      while (loglik_new < loglik && step > 1e-10)
      {
        step = step / 2;
        coef_new = coef + step * desend_direction;
        i = 1;
        while(i < k)
        {
          for (i = 1; i < k; i++)
          {
            if (coef_new(i) <= coef_new(i - 1))
            {
              step = step / 2;
              coef_new = coef + step * desend_direction;
              break;
            }
          }
        }
        beta.col(0) = coef_new.tail(p);
        coef0.head(k) = coef_new.head(k);
        loglik_new = -loss_function(X, y, weights, beta, coef0, lambda);
      }


      printf("%d: g = %f, d_d = %f, step = %.10f, loss = %f\n", j, g.cwiseAbs2().sum(), desend_direction.cwiseAbs2().sum(), step, -loglik_new);

      coef = coef_new;
      loglik = loglik_new;
    }

    return coef;
  }



