#ifndef QUERY_H
#define QUERY_H

#if defined __cplusplus

#include "../inst/include/iscream_types.h"
#include <Rcpp.h>
#include <string>
#include "parsers.hpp"

std::vector<std::string> tabix_query(const std::string& region, htsFile* bedFile, tbx_t* tbx);
std::vector<std::vector<std::string>> query_intervals(const char* fname, std::vector<std::string>& regions);
Rcpp::List query_regions_from_file(const char* fname, std::vector<std::string>& regions);
std::vector<std::vector<std::string>> query_interval(std::vector<std::string>& bedfiles, std::string& region);
std::vector<std::string> query_interval(std::string& bedfile, std::string& region);

#endif /* __cplusplus */

#endif /* ifndef QUERY_H */
