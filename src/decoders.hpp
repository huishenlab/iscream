#ifndef DECODERS_H
#define DECODERS_H

#if defined __cplusplus

#include <Rcpp.h>

int decode_beta(uint32_t encoded);
int decode_cov(uint32_t encoded);
int decode_m(uint32_t encoded);
Rcpp::IntegerVector vdecoder(Rcpp::IntegerVector& encoded, int measure);
Rcpp::DoubleVector vdouble_decoder(Rcpp::DoubleVector& encoded, int measure);

#endif /* __cplusplus */

#endif /* ifndef DECODERS_H */
