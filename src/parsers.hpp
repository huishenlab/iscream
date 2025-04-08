#ifndef PARSERS_H
#define PARSERS_H

# if defined __cplusplus

#include <string>
#include <Rcpp.h>

struct BedLine {
    std::string chr;
    int start;
    int end;
    float beta;
    int cov;
    int m_count;
};

std::vector<std::string_view> split_bedstring(std::string_view bedString);
BedLine parseBiscuitRecord(const std::string& bedString);
BedLine parseCovRecord(const std::string& bedString);
/*std::string CpGID(BedLine& parsed_bedline);*/

#endif /* __cplusplus */

#endif /* ifndef PARSERS_H */
