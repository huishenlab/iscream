#ifndef PARSERS_H
#define PARSERS_H

# if defined __cplusplus

#include <Rcpp.h>
#include <string>
#include "encoders.hpp"

struct BedLine {
    std::string chr;
    int start;
    int end;
    int beta;
    int cov;
};

};

BedLine parseBEDRecord(const std::string& bedString);

#endif /* __cplusplus */

#endif /* ifndef PARSERS_H */
