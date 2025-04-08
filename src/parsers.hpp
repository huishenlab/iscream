#ifndef PARSERS_H
#define PARSERS_H

# if defined __cplusplus

#include <string>
#include <Rcpp.h>

struct BedRecord {
    std::string chr;
    int start;
    int end;
    float data[3];
};

std::vector<std::string_view> split_bedstring(std::string_view bedString);
BedLine parseBiscuitRecord(const std::string& bedString);
BedLine parseCovRecord(const std::string& bedString);
/*std::string CpGID(BedLine& parsed_bedline);*/

#endif /* __cplusplus */

#endif /* ifndef PARSERS_H */
