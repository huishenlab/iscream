#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
int encode_meth(uint32_t meth, uint32_t unmeth) {
  return (meth << 16) + unmeth;
}

// [[Rcpp::export]]
int decode_meth(uint32_t encoded) {
  return encoded >> 16;
}

// [[Rcpp::export]]
int decode_unmeth(uint32_t encoded) {
  return (encoded) & 0xffff;
}

// [[Rcpp::export]]
int encoder(float beta, uint32_t cov) {
    int encoded = std::round(beta * 100) + (cov << 16);
    return encoded;
}

// [[Rcpp::export]]
IntegerVector decode_beta(IntegerVector& beta_col) {
    for (int i = 0; i < beta_col.length(); i++) {
        beta_col[i] = (beta_col[i]) & 0xffff;
    }
    return beta_col;
}
