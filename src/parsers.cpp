#include "parsers.hpp"

//' Split a line from a bed file
//'
//' @param A line from a bed file
std::vector<std::string_view> split_bedstring(std::string_view bedString) {
    // https://www.reddit.com/r/Cplusplus/comments/gnc9rz/fast_string_split/
    std::vector<std::string_view> res;
    res.reserve(bedString.length() / 2);

    const char* ptr = bedString.data();
    size_t size = 0;

    for(const char c : bedString) {
        if(c == '\t') {
            res.emplace_back(ptr, size);
            ptr += size + 1;
            size = 0;
        }
        else {
            ++size;
        }
    }

    if(size)
        res.emplace_back(ptr, size);
    return res;
}

//' Parse a bed record into chr, start, and compressed beta and cov
//'
//' @param A line from a bed file
BedLine parseBEDRecord(const std::string& bedString) {
    std::vector<std::string_view> res = split_bedstring(bedString);
    std::vector<std::string> fields(res.begin(), res.end());
    std::string chrom = fields[0];
    int start = std::stoi(fields[1]);
    int end = std::stoi(fields[2]);
    float beta = std::stof(fields[3]);
    int cov = std::stoi(fields[4]);
    int m_count = (int) std::round(cov * beta);

    return BedLine{chrom, start, end, beta, cov, m_count};
}

//' Parse a bed record into chr, start, and compressed beta and cov
//'
//' @param A line from a bed file
BedLine parseCovRecord(const std::string& bedString) {
    std::vector<std::string_view> res = split_bedstring(bedString);
    // TODO: unmapped memory if there are fewer than 5 detected fields
    std::vector<std::string> fields(res.begin(), res.end());
    std::string chrom = fields[0];
    int start = std::stoi(fields[1]);
    int end = std::stoi(fields[2]);
    int m_count = std::stoi(fields[4]);
    int u_count = std::stoi(fields[5]);
    int cov = m_count + u_count;
    float beta = m_count / (m_count + u_count);

    BedLine read = {chrom, start, end, beta, cov, m_count};
    return read;
}

EncodedBedLine encodeBedRecord(const std::string& bedString) {
    BedLine parsed_bedline = parseBEDRecord(bedString);
    int encoded = encoder(parsed_bedline.beta, parsed_bedline.cov);
    return EncodedBedLine{parsed_bedline.chr, parsed_bedline.start, encoded};
}

std::string CpGID(BedLine& parsed_bedline) {

    std::stringstream cpg_stream;

    cpg_stream << parsed_bedline.chr << ":" << parsed_bedline.start;
    return cpg_stream.str();
}
