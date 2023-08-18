#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export]]
int decode_beta(uint32_t encoded) {
    return (int) encoded & 0xffff;
}

//' @export
// [[Rcpp::export]]
int decode_cov(uint32_t encoded) {
    return (int) encoded >> 16;
}

//' Vector Decoder
//' This function unpacks the encoded \eqn{\beta} and coverage values from the
//' encoded column of a methylation matrix.
//' @param encoded An encoded sample column from a data.table methylation matrix of CpG loci by samples
//' @param measure The required decoded value: 1 for beta, 2 for coverage and 3 for M
//' @export
// [[Rcpp::export]]
IntegerVector vdecoder(IntegerVector& encoded, int measure) {
  for (int i = 0; i < encoded.length(); i++) {
    if (IntegerVector::is_na(encoded[i])) {
      encoded[i] = NA_INTEGER;
    } else if (measure == 1) {
        encoded[i] = decode_beta(encoded[i]);
    } else if (measure == 2) {
        encoded[i] = decode_cov(encoded[i]);
    } else if (measure == 3) {
        encoded[i] = std::round((float) (decode_beta(encoded[i]) * decode_cov(encoded[i])) / 100);
    }
  }
  return encoded;
}


//' Double Vector Decoder
//' This function unpacks the encoded \eqn{\beta} and coverage values from the
//' encoded column of a methylation matrix.
//' @param encoded An encoded sample column from a data.table methylation matrix of CpG loci by samples
//' @param measure The required decoded value: 1 for beta, 2 for coverage and 3 for M
//' @export
// [[Rcpp::export]]
DoubleVector vdouble_decoder(DoubleVector& encoded, int measure) {
  for (int i = 0; i < encoded.length(); i++) {
    if (DoubleVector::is_na(encoded[i])) {
      encoded[i] = NA_INTEGER;
    } else if (measure == 1) {
        encoded[i] = decode_beta(encoded[i]);
    } else if (measure == 2) {
        encoded[i] = decode_cov(encoded[i]);
    } else if (measure == 3) {
        encoded[i] = std::round((float) (decode_beta(encoded[i]) * decode_cov(encoded[i])) / 100);
    }
  }
  return encoded;
}

