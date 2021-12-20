#include <Rcpp.h>
using namespace Rcpp;


//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp, the target density has been given.
//' @param num The number of samples. Default is \code{10000}.
//' @param a Parameter of density. Default is \code{1}.
//' @param b Parameter of density. Default is \code{1}.
//' @param n Parameter of density. Default is \code{100}.
//' @return A \code{matrix}, the number of rows is the param num and the numberof colmuns is 2.
//' @examples
//' \dontrun{
//' c.chain <- bi_chain_cpp()
//' plot(c.chain)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix bi_chain_cpp(int num=10000,int a=1,int b=1,int n=100){
    NumericMatrix chain(num,2);
    chain(0,0) = floor(n/2);
    chain(0,1) = 0.5;

    for (int i=1; i<num; i++) {
        chain(i, 0) = rbinom(1,n,chain(i-1, 1))[0];
        chain(i, 1) = rbeta(1,(int)chain(i, 0)+a,n-(int)chain(i, 0)+b)[0];
    }

    return chain;
}



