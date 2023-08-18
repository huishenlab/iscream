#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

static inline int encoder(float beta, int cov) {
    return ((int) std::round(beta * 100) + (cov << 16));
}

//' Vector Encoder
//' Bit-packs a vector of beta and coverage values
//' @param beta_col Name of column containing beta values
//' @param cov_col Name of column containing coverage values
//' @export
// [[Rcpp::export]]
IntegerVector vencoder(NumericVector beta_col, IntegerVector cov_col) {
  IntegerVector encoded_col (beta_col.length());

  for (int i = 0; i < beta_col.length(); i++) {
    encoded_col[i] = encoder(beta_col[i], cov_col[i]);
  }
  return encoded_col;
}

