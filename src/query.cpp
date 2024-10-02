#include "query.hpp"

//' Query a genomic interval from a opened htsFile and return the reads in it
//'
//' @param region Genomic region string in the form "chr:start-end"
//' @param bedFile The opened htslib bed file stream
//' @param tbx The bed file's tab-index
//' @returns A vector of strings from the matching region of the bed file
std::vector<std::string> tabix_query(
    const std::string& region,
    htsFile* bedFile, tbx_t* tbx
) {
    std::vector<std::string> reads;
    hts_itr_t* iter = tbx_itr_querys(tbx, region.c_str());
    // TODO: catch errors here
    // its currently just returning nothing if the region string is invalid
    // Rcpp::stop causes a core dump instead of exiting gracefully because of
    // openmp threads dying before being destroyed
    kstring_t str = {0, 0, 0};

    while (tbx_itr_next(bedFile, tbx, iter, &str) >= 0) {
        reads.push_back(str.s);
    }

    free(str.s);
    hts_itr_destroy(iter);

    return reads;
}

//' Get reads from multiple genomic regions from a tabixed bed file
//'
//' @param bedfile The name of the bed file - must have a corresponding tabix
//' file with the same name and .tbi extension
//' @param regions A vector of region strings in the form "chr:start-end"
std::vector<RegionQuery> query_intervals(
    const char* bedfile,
    const std::vector<std::string>& regions
) {

    htsFile* bedFile = hts_open(bedfile, "r");
    tbx_t* tbx = tbx_index_load3(bedfile, NULL, 0);

    std::vector<RegionQuery> all_reads(regions.size());

    for (int i = 0; i < regions.size(); i++) {
        all_reads[i] = RegionQuery(regions[i], tabix_query(regions[i], bedFile, tbx));
    }

    tbx_destroy(tbx);
    hts_close(bedFile);

    return all_reads;
}

//' Get reads from single genomic regions from multiple tabixed bed file.
//'
//' @param bedfiles A vector of bedfile names - must have corresponding tabix
//' files with the same name and .tbi extension
//' @param region A vector regions string in the form "chr:start-end"
std::vector<std::vector<std::string>> query_interval(
    const std::vector<std::string>& bedfiles,
    const std::string& region
) {

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

//' Get reads from a single genomic region from one tabixed bed file.
//'
//' @param bedfile The name of the bed file - must have a corresponding tabix
//' file with the same name and .tbi extension
//' @param region The region string in the form "chr:start-end"
std::vector<std::string> query_interval(
    const std::string& bedfile,
    const std::string& region
) {
    htsFile* bedFile = hts_open(bedfile.c_str(), "r");
    tbx_t* tbx = tbx_index_load3(bedfile.c_str(), NULL, 0);

    std::vector<std::string> regional_cpgs = tabix_query(region, bedFile, tbx);

    tbx_destroy(tbx);
    hts_close(bedFile);

    return regional_cpgs;
}

//' Query the chromosomes or seqnames from a vector of files
//' @param bedfile_vec The vector of bedfile paths
//' @return A vector of seqnames
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
std::set<std::string> Cpp_query_chroms(const std::vector<std::string>& bedfile_vec, const int nthreads) {
    std::set<std::string> seqnames;

    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for (const std::string fname : bedfile_vec) {
        tbx_t* tbx = tbx_index_load3(fname.c_str(), NULL, 0);
        if (!tbx) spdlog::error("Could not load .tbi index of {}", fname);
        int i, nseq = 1<<1;
        const char **seq = tbx_seqnames(tbx, &nseq);
        if (!seq) spdlog::error("Could not get list of chromosome names");
        for (int i = 0; i < nseq; i++) {
            #pragma omp critical
            {
                seqnames.insert(seq[i]);
            }
        }
    }
    return seqnames;
}

