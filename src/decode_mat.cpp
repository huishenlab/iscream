// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Beta value decoder
//'
//' @param encoded The bit-packed beta and cov Int
//' @return The beta value
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
int decode_beta(int encoded) {
    return (int) encoded & 0xffff;
}

//' Coverage value decoder
//'
//' @param encoded The bit-packed beta and cov Int
//' @return The coverage value
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
int decode_cov(int encoded) {
    return (int) encoded >> 16;
}

//' @keywords internal
//' @export
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

//' @keywords internal
//' @export
// [[Rcpp::export]]
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


//' @keywords internal
//' @export
// [[Rcpp::export]]
void get_beta(Rcpp::NumericMatrix& m, const int nthreads) {
  arma::mat M(m.begin(), m.nrow(), m.ncol(), false);

     #if defined(_OPENMP)
         #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int coln = 0; coln < M.n_cols; coln++) {
        M.col(coln).transform( [](double packed_val) { return ( (int) packed_val & 0xffff); } );
    }
}
