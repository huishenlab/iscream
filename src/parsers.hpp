#ifndef PARSERS_H
#define PARSERS_H

# if defined __cplusplus

#include <Rcpp.h>
#include <string>
#include "encoders.hpp"

struct BedLine {
    const char* chr;
    int start;
    int encoded;

};

BedLine parseBEDRecord(const std::string bedString);

#endif /* __cplusplus */

#endif /* ifndef PARSERS_H */
