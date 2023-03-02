#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int encode_meth(uint32_t meth, uint32_t unmeth) {
  return (meth << 16) + unmeth;
}

//' @export
// [[Rcpp::export]]
int decode_meth(uint32_t encoded) {
  return encoded >> 16;
}

//' @export
// [[Rcpp::export]]
int decode_unmeth(uint32_t encoded) {
  return (encoded) & 0xffff;
}

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

