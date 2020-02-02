#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat libNormFast(arma::mat raw_matrix, arma::rowvec scalefactor)
  {
    arma::rowvec libsizes = arma::sum(raw_matrix, 0);
    raw_matrix.each_row() %= 1/libsizes;
    raw_matrix.each_row() %= scalefactor;
    return raw_matrix;
  }

// [[Rcpp::export]]
arma::mat logNormFast(arma::mat mymatrix)
  {
    arma::mat lognormmatrix = arma::log1p(mymatrix);
    return lognormmatrix;
  }

// [[Rcpp::export]]
arma::mat armaScaleRow(arma::mat Z)
  {
  arma::mat Zt = Z.t();
  unsigned int j, n = Zt.n_rows, k = Zt.n_cols;
  double avg, sd;
  arma::colvec z;
  arma::mat res = arma::zeros(n, k);

  for (j=0; j<k; j++) {
    z = Zt.col(j);
    avg = arma::mean(z);
    sd = arma::stddev(z);
    res.col(j) = (z - avg) / sd;
  }

  return res.t();
}

// [[Rcpp::export]]
arma::mat armaScale(arma::mat Z)
  {
  unsigned int j, n = Z.n_rows, k = Z.n_cols;
  double avg, sd;
  arma::colvec z;
  arma::mat res = arma::zeros(n, k);

  for (j=0; j<k; j++) {
    z = Z.col(j);
    avg = arma::mean(z);
    sd = arma::stddev(z);
    res.col(j) = (z - avg) / sd;
  }

  return res;
}





