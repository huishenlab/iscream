#ifndef QUERY_H
#define QUERY_H

#if defined __cplusplus

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#include <htslib/hts.h>
#include <htslib/tbx.h>
#include <string>
#include "parsers.hpp"
#include "log.hpp"
#include "set"


struct RegionQuery{
    RegionQuery() : interval_str(), cpgs_in_interval() {}
    RegionQuery(std::string interval_str, std::vector<std::string> cpgs_in_interval) : interval_str(interval_str), cpgs_in_interval(cpgs_in_interval){}
    std::string interval_str;
    std::vector<std::string> cpgs_in_interval;
};

typedef std::vector<RegionQuery> MultiRegionQuery;


std::vector<std::string> tabix_query(const std::string& region, htsFile* bedFile, tbx_t* tbx);
std::vector<RegionQuery> query_intervals(const char* fname, const std::vector<std::string>& regions);
std::vector<std::vector<std::string>> query_interval(const std::vector<std::string>& bedfiles, const std::string& region);
std::vector<std::string> query_interval(const std::string& bedfile, const std::string& region);
std::set<std::string> Cpp_query_chroms(const std::vector<std::string>& bedfile_vec, const int nthreads);

#endif /* __cplusplus */

#endif /* ifndef QUERY_H */
