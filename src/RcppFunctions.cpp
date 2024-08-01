#include "bsseq.hpp"

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
Rcpp::List Cpp_query_all(std::vector<std::string>& bedfiles, std::vector<std::string>& regions, const bool bismark, const bool merged, const int nthreads) {

    BS bsseq = BS(bedfiles, regions, bismark, merged, nthreads);
    return bsseq.wrap();
}
