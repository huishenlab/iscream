#include "bsseq.hpp"

//' Make a Bsseq object
//' @param bedfiles A vector of bedfiles
//' @param regions A vector of regions
//' @param print Whether to print the Cov and M matrices
//' @export
// [[Rcpp::export]]
Rcpp::S4 make_bsseq(std::vector<std::string>& bedfiles, std::vector<std::string>& regions, bool print) {

    BS bsseq = BS(bedfiles, regions);
    if (print) {
        bsseq.print_BS();
    }
    return bsseq.wrap();
}
