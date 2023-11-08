#include "query.hpp"

//' Query a genomic interval from a opened htsFile and return the reads in it
//' @param region Genomic region string in the form "chr:start-end"
//' @param bedFile The opened htslib bed file stream
//' @param tbx The bed file's tab-index
//' @returns A vector of strings from the matching region of the bed file
std::vector<std::string> tabix_query(const std::string& region, htsFile* bedFile, tbx_t* tbx) {
    std::vector<std::string> reads;
    hts_itr_t* iter = tbx_itr_querys(tbx, region.c_str());
    kstring_t str = {0, 0, 0};

    while (tbx_itr_next(bedFile, tbx, iter, &str) >= 0) {
        reads.push_back(str.s);
    }

    free(str.s);
    hts_itr_destroy(iter);

    return reads;
}

//' Get reads from multiple genomic regions from a tabixed bed file
//' @param fname The name of the bed file - must have a tabix file with the same name and .tbi extension
//' @param regions A vector of regions strings of the form "chr:start-end"
std::vector<std::vector<std::string>> query_intervals(const char* fname, std::vector<std::string>& regions) {

    htsFile* bedFile = hts_open(fname, "r");
    tbx_t* tbx = tbx_index_load3(fname, NULL, 0);

    std::vector<std::vector<std::string>> all_reads(regions.size());

    for (int i = 0; i < regions.size(); i++) {
        all_reads[i] = tabix_query(regions[i], bedFile, tbx);
    }

    tbx_destroy(tbx);
    hts_close(bedFile);

    return all_reads;
}

//' Get list of reads from multiple genomic regions from a tabixed bed file.
//' @param fname The name of the bed file - must have a tabix file with the same name and .tbi extension
//' @param regions A vector of regions strings of the form "chr:start-end"
//' @export
// [[Rcpp::export]]
Rcpp::List query_regions_from_file(const char* fname, std::vector<std::string>& regions) {

    std::vector<std::vector<std::string>> reads_in_file = query_intervals(fname, regions);
    Rcpp::List read_list = Rcpp::wrap(reads_in_file);

    read_list.attr("names") = regions;

    return read_list;
}

//' Get reads from single genomic regions from multiple tabixed bed file.
//' @param fname The name of the bed file - must have a tabix file with the same name and .tbi extension
//' @param regions A vector of regions strings of the form "chr:start-end"
//' @export
// [[Rcpp::export]]
std::vector<std::vector<std::string>> query_interval(std::vector<std::string>& bedfiles, std::string& region) {

    std::vector<std::vector<std::string>> all_reads(bedfiles.size());

    for (int i = 0; i < bedfiles.size(); i++) {
        htsFile* bedFile = hts_open(bedfiles[i].c_str(), "r");
        tbx_t* tbx = tbx_index_load3(bedfiles[i].c_str(), NULL, 0);

        all_reads.push_back(tabix_query(region, bedFile, tbx));

        tbx_destroy(tbx);
        hts_close(bedFile);
    }

    return all_reads;

}
