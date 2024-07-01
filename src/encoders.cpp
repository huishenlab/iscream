#include "encoders.hpp"

using namespace Rcpp;

//' Vector Encoder
//'
//' Bit-packs a vector of beta and coverage values
//' @param beta_col Name of column containing beta values
//' @param cov_col Name of column containing coverage values
//'
//' @keywords experimental
//' @export
// [[Rcpp::export]]
IntegerVector vencoder(NumericVector beta_col, IntegerVector cov_col) {
  IntegerVector encoded_col (beta_col.length());

  for (int i = 0; i < beta_col.length(); i++) {
    encoded_col[i] = encoder(beta_col[i], cov_col[i]);
  }
  return encoded_col;
}

