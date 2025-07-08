// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Beta value decoder
//
// @param encoded The bit-packed beta and cov Int
// @return The beta value
//
// [[Rcpp::export]]
int decode_beta(int encoded) {
    return (int) encoded & 0xffff;
}

// Coverage value decoder
//
// @param encoded The bit-packed beta and cov Int
// @return The coverage value
//
// [[Rcpp::export]]
int decode_cov(int encoded) {
    return (int) encoded >> 16;
}

// Get coverage values from bit-packed dense matrix
//
// [[Rcpp::export]]
void get_cov(Rcpp::NumericMatrix& m, const int nthreads) {
  arma::mat M(m.begin(), m.nrow(), m.ncol(), false);

     #if defined(_OPENMP)
         #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int coln = 0; coln < M.n_cols; coln++) {
        M.col(coln).transform( [](double packed_val) { return ( (int) packed_val >> 16); } );
    }
}

// Get M values from bit-packed dense matrix
//
//[[Rcpp::export]]
void get_m(Rcpp::NumericMatrix& m, const int nthreads) {
  arma::mat M(m.begin(), m.nrow(), m.ncol(), false);

     #if defined(_OPENMP)
         #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int coln = 0; coln < M.n_cols; coln++) {
        M.col(coln).transform( [](double packed_val) {
            return std::round((float) (decode_beta(packed_val) * decode_cov(packed_val)) / 100);
        } );
    }
}


// Get beta values from bit-packed dense matrix
//
// [[Rcpp::export]]
void get_beta(Rcpp::NumericMatrix& m, const int nthreads) {
  arma::mat M(m.begin(), m.nrow(), m.ncol(), false);

     #if defined(_OPENMP)
         #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int coln = 0; coln < M.n_cols; coln++) {
        M.col(coln).transform( [](double packed_val) { return ( (double) ((int) packed_val & 0xffff) / 100); } );
    }
}

// Get coverage values from bit-packed sparse matrix
//
// [[Rcpp::export]]
void get_cov_sparse(Rcpp::S4& m) {
    Rcpp::NumericVector elements = m.slot("x");
    arma::vec x(elements.begin(), elements.size(), false);

    for (int i = 0; i < x.size(); i++) {
        x[i] = (int) x[i] >> 16;
    }
}

// Get M values from bit-packed sparse matrix
//
// [[Rcpp::export]]
void get_m_sparse(Rcpp::S4& m) {
    Rcpp::NumericVector elements = m.slot("x");
    arma::vec x(elements.begin(), elements.size(), false);

    for (int i = 0; i < x.size(); i++) {
        x[i] = std::round((float) (decode_beta(x[i]) * decode_cov(x[i])) / 100);
    }
}


// Get beta values from bit-packed sparse matrix
//
// @param m The bit-packed dgcMatrix
// @return The beta value
//
// [[Rcpp::export]]
void get_beta_sparse(Rcpp::S4& m) {
    Rcpp::NumericVector elements = m.slot("x");
    arma::vec x(elements.begin(), elements.size(), false);

    for (int i = 0; i < x.size(); i++) {
        x[i] = (double) ((int) x[i] & 0xffff) / 100;
    }
}
