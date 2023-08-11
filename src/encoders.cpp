#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

static inline int encoder(float beta, uint32_t cov) {
    return (int)(std::round(beta * 100) + (cov << 16));
}

//' @export
// [[Rcpp::export]]
IntegerVector vencoder(NumericVector beta_col, IntegerVector cov_col) {
  IntegerVector encoded_col (beta_col.length());

  for (int i = 0; i < beta_col.length(); i++) {
    encoded_col[i] = encoder(beta_col[i], cov_col[i]);
  }
  return encoded_col;
}

//' @export
// [[Rcpp::export]]
IntegerVector decode_beta(IntegerVector& beta_col) {
    for (int i = 0; i < beta_col.length(); i++) {
        beta_col[i] = (beta_col[i]) & 0xffff;
    }
    return beta_col;
}

