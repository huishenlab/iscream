#include "bsseq.hpp"

//' Make a Bsseq object
//' @param bedfiles A vector of bedfiles
//' @param regions A vector of regions
//' @export
// [[Rcpp::export]]
Rcpp::S4 make_bsseq(std::vector<std::string>& bedfiles, std::vector<std::string>& regions) {

    BS bsseq = BS(bedfiles, regions);
    return bsseq.wrap();
}
