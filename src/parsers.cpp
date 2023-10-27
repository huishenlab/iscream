#include "parsers.hpp"

//' Parse a bed record into chr, start, and compressed beta and cov
//' @param A line from a bed file
BedLine parseBEDRecord(const std::string& bedString) {
    std::istringstream ss(bedString);
    std::string token;
    std::vector<std::string> fields;

    while (std::getline(ss, token, '\t')) {
        fields.push_back(token);
    }

    std::string chrom = fields[0];
    int start = std::stoi(fields[1]);
    int encoded = encoder(std::stof(fields[3]), std::stoi(fields[4]));

    BedLine read = {chrom.c_str(), start, encoded};

    return read;
}

