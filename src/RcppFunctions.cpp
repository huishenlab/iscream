#include "bsseq.hpp"

//' Make a Bsseq object
//' @param bedfiles A vector of bedfiles
//' @param regions A vector of regions
//' @param nthreads Set number of threads to use overriding the
//' `"iscream.threads"` option. See `?set_threads` for more information.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::S4 Cpp_query_all(std::vector<std::string>& bedfiles, std::vector<std::string>& regions, const int nthreads) {

    BS bsseq = BS(bedfiles, regions, nthreads);
    return bsseq.wrap();
}
