#include "query_all.hpp"

//' Make the inputs to a BSseq object
//' Produces M and coverage matrices, 'chr:start' of the loci as rownames, sample names
//' @param bedfiles A vector of bedfiles
//' @param regions A vector of regions
//' @param bismark Whether the input is a bismark coverage file
//' @param prealloc The number of rows to initialize the matrices with
//' @param nthreads Set number of threads to use overriding the
//' `"iscream.threads"` option. See `?set_threads` for more information.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List Cpp_query_all_dense(
    std::vector<std::string>& bedfiles,
    std::vector<std::string>& regions,
    const bool bismark,
    const bool merged,
    const bool sparse,
    const int prealloc,
    const int nthreads) {

    QueryAll bsseq = QueryAll<arma::umat>(bedfiles, regions, bismark, merged, sparse, prealloc, nthreads);
    return bsseq.wrap();
}

//' Make the inputs to a BSseq object
//' Produces M and coverage matrices, 'chr:start' of the loci as rownames, sample names
//' @param bedfiles A vector of bedfiles
//' @param regions A vector of regions
//' @param bismark Whether the input is a bismark coverage file
//' @param nthreads Set number of threads to use overriding the
//' `"iscream.threads"` option. See `?set_threads` for more information.
//'
//' @export
// [[Rcpp::export]]
Rcpp::List Cpp_query_all_sparse(
    std::vector<std::string>& bedfiles,
    std::vector<std::string>& regions,
    const bool bismark,
    const bool merged,
    const bool sparse,
    const int prealloc,
    const int nthreads) {

    QueryAll bsseq = QueryAll<arma::sp_umat>(bedfiles, regions, bismark, merged, sparse, prealloc, nthreads);
    return bsseq.wrap();
}

//' Make the inputs to a BSseq object
//' Produces M and coverage matrices, 'chr:start' of the loci as rownames,
//' sample names
//' @param bedfiles A vector of bedfiles
//' @param mmat A NumericMatrix from R for M values - must be large enough to
//' hold the final matrix produced or a copy will be made
//' @param cmat A NumericMatrix from R for coverage - must be large enough to
//' hold the final matrix produced or a copy will be made
//' @param regions A vector of regions
//' @param bismark Whether the input is a bismark coverage file
//' @param nthreads Set number of threads to use overriding the
//' `"iscream.threads"` option. See `?set_threads` for more information.
//'
//' @export
// [[Rcpp::export]]
void Cpp_query_all_reference(
    std::vector<std::string>& bedfiles,
    std::vector<std::string>& regions,
    Rcpp::NumericMatrix& mmat,
    Rcpp::NumericMatrix& cmat,
    const bool bismark,
    const bool merged,
    const bool sparse,
    const int nthreads) {

    QueryAll bsseq = QueryAll<arma::mat>(bedfiles, regions, mmat, cmat, bismark, merged, sparse, nthreads);
}
