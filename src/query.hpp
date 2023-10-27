#ifndef QUERY_H
#define QUERY_H

#if defined __cplusplus

#include "../inst/include/iscream_types.h"
#include <Rcpp.h>
#include <string>
#include "parsers.hpp"

std::vector<std::string> query_region(const std::string& region, htsFile* bedFile, tbx_t* tbx);
#endif /* __cplusplus */

#endif /* ifndef QUERY_H */
