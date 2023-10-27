#ifndef QUERY_H
#define QUERY_H

#if defined __cplusplus

#include "../inst/include/iscream_types.h"
#include <Rcpp.h>
#include <string>
#include "parsers.hpp"

std::vector<std::string> query_region(const std::string& region, htsFile* bedFile, tbx_t* tbx);
std::vector<std::vector<std::string>> query_file_priv(const char* fname, std::vector<std::string>& regions);
Rcpp::List query_file(const char* fname, std::vector<std::string>& regions);

#endif /* __cplusplus */

#endif /* ifndef QUERY_H */
