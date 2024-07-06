#ifndef PARSERS_H
#define PARSERS_H

# if defined __cplusplus

#include <string>
#include "encoders.hpp"

struct BedLine {
    std::string chr;
    int start;
    int end;
    float beta;
    int cov;
    int m_count;
};

struct EncodedBedLine {
    std::string chr;
    int start;
    int encoded;
};

BedLine parseBEDRecord(const std::string& bedString);
BedLine parseCovRecord(const std::string& bedString);
EncodedBedLine encodeBedRecord(const std::string& bedString);
std::string CpGID(BedLine& parsed_bedline);

#endif /* __cplusplus */

#endif /* ifndef PARSERS_H */
