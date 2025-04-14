#ifndef PARSERS_H
#define PARSERS_H

# if defined __cplusplus

#include <string>
#include <Rcpp.h>

struct BedRecord {
    std::string chr;
    int start;
    int end;
    std::vector<float> data;
    size_t size;
};

std::vector<std::string_view> split_bedstring(std::string_view bedString);
/*std::string CpGID(BedLine& parsed_bedline);*/
BedRecord parseBedRecord(const std::string& bedString, const int valInd1);
BedRecord parseBedRecord(const std::string& bedString, std::vector<int> valInd);
BedRecord parseBiscuitRecord(const std::string& bedString);
BedRecord parseCovRecord(const std::string& bedString);

#endif /* __cplusplus */

#endif /* ifndef PARSERS_H */
