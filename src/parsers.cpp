#include "parsers.hpp"

//' Parse a bed record into chr, start, and compressed beta and cov
//'
//' @param A line from a bed file
BedLine parseBEDRecord(const std::string& bedString) {
    std::istringstream ss(bedString);
    std::string token;
    std::vector<std::string> fields(5);
    int index = 0;

    while (std::getline(ss, token, '\t') && index < fields.size()) {
        fields[index] = token;
        index++;
    }

    // TODO: unmapped memory if there are fewer than 5 detected fields
    std::string chrom = fields[0];
    int start = std::stoi(fields[1]);
    int end = std::stoi(fields[2]);
    float beta = std::stof(fields[3]);
    int cov = std::stoi(fields[4]);
    int m_count = (int) std::round(cov * beta);

    BedLine read = {chrom, start, end, beta, cov, m_count};
    return read;
}

//' Parse a bed record into chr, start, and compressed beta and cov
//'
//' @param A line from a bed file
BedLine parseCovRecord(const std::string& bedString) {
    std::istringstream ss(bedString);
    std::string token;
    std::vector<std::string> fields(6);
    int index = 0;

    while (std::getline(ss, token, '\t') && index < fields.size()) {
        fields[index] = token;
        index++;
    }

    // TODO: unmapped memory if there are fewer than 5 detected fields
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
