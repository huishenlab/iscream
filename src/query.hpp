#ifndef QUERY_H
#define QUERY_H

#if defined __cplusplus

#include "../inst/include/iscream_types.h"
#include <string>
#include "parsers.hpp"

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
std::vector<std::string> query_chroms(const std::string& fname);

#endif /* __cplusplus */

#endif /* ifndef QUERY_H */
