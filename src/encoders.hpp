#ifndef ENCODERS_H
#define ENCODERS_H

#if defined __cplusplus

#include <Rcpp.h>

inline int encoder(float beta, int cov) {
    return ((int) std::round(beta * 100) + (cov << 16));
}

Rcpp::IntegerVector vencoder(Rcpp::NumericVector beta_col, Rcpp::IntegerVector cov_col);

#endif /* __cplusplus */

#endif /* ifndef ENCODERS_H */
