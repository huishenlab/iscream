#include "query_all.hpp"

//' Make a Bsseq object
//' @param bedfiles A vector of bedfiles
//' @param regions A vector of regions
//' @param bismark Whether the input is a bismark coverage file
//' @param nthreads Set number of threads to use overriding the
//' `"iscream.threads"` option. See `?set_threads` for more information.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List Cpp_query_all(std::vector<std::string>& bedfiles, std::vector<std::string>& regions, const bool bismark, const bool merged, const bool sparse, const int nthreads) {

    if (sparse) {
        QueryAll bsseq = QueryAll<arma::sp_umat>(bedfiles, regions, bismark, merged, sparse, nthreads);
        return bsseq.wrap();
    } else {
        QueryAll bsseq = QueryAll<arma::umat>(bedfiles, regions, bismark, merged, sparse, nthreads);
        return bsseq.wrap();
    }
}
