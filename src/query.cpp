#include "query.hpp"

//' Query a genomic interval and return the reads in it
//' @param region Genomic region string in the form "chr:start-end"
//' @param bedFile The opened htslib bed file stream
//' @param tbx The bed file's tab-index
//' @returns A vector of strings from the matching region of the bed file
std::vector<std::string> query_region(const std::string& region, htsFile* bedFile, tbx_t* tbx) {
    std::vector<std::string> reads;
    hts_itr_t* iter = tbx_itr_querys(tbx, region.c_str());
    kstring_t str = {0, 0, 0};

    while (tbx_itr_next(bedFile, tbx, iter, &str) >= 0) {
        reads.push_back(str.s);
        /* printf("%s\n", str.s); */
    }

    free(str.s);
    hts_itr_destroy(iter);

    return reads;
}

