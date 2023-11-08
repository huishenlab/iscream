#include <cmath>
#include <Rcpp.h>
#include "decoders.hpp"

using namespace Rcpp;

//' Beta value decoder
//' @param encoded The bit-packed beta and cov Int
//' @return The beta value
//' @export
// [[Rcpp::export]]
int decode_beta(uint32_t encoded) {
    return (int) encoded & 0xffff;
}

//' Coverage value decoder
//' @param encoded The bit-packed beta and cov Int
//' @return The coverage value
//' @export
// [[Rcpp::export]]
int decode_cov(uint32_t encoded) {
    return (int) encoded >> 16;
}

//' M value decoder
//' @param encoded The bit-packed beta and cov Int
//' @return The M value
//' @export
// [[Rcpp::export]]
int decode_m(uint32_t encoded) {
    return std::round((float) (decode_beta(encoded) * decode_cov(encoded)) / 100);
}

//' Vector Decoder
//' This function unpacks the encoded \eqn{\beta} and coverage values from the
//' encoded column of a methylation matrix.
//' @param encoded An encoded sample column from a data.table methylation matrix of CpG loci by samples
//' @param measure The required decoded value: 1 for beta, 2 for coverage and 3 for M
//' @return An IntegerVector of unpacked beta, coverage or M values
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
        encoded[i] = decode_m(encoded[i]);
    }
  }
  return encoded;
}


//' Double Vector Decoder
//' This function unpacks the encoded \eqn{\beta} and coverage values from the
//' encoded column of a methylation matrix. Same as `vdecoder`, but accepts double columns
//' @param encoded An encoded sample column from a data.table methylation matrix of CpG loci by samples
//' @param measure The required decoded value: 1 for beta, 2 for coverage and 3 for M
//' @return A DoubleVector of unpacked beta, coverage or M values
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
        encoded[i] = decode_m(encoded[i]);
    }
  }
  return encoded;
}

