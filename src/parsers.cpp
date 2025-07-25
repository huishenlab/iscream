#include "parsers.hpp"
#include "log.hpp"

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

BedRecord parseBedRecord(const std::string& bedString, std::vector<int> valInd) {
    std::vector<std::string_view> res = split_bedstring(bedString);
    std::vector<std::string> fields(res.begin(), res.end());
    size_t col_count = valInd.size();

    if (col_count >= fields.size()) {
        throw std::out_of_range(fmt::format("Error: Column {} does not exist\n", col_count));
    }

    std::string chrom = fields[0];
    int start = std::stoi(fields[1]);
    int end = std::stoi(fields[2]);
    std::vector<float> data(col_count);
    int i;
    try {
        for (i = 0; i < valInd.size(); i++) {
            data[i] = std::stof(fields[valInd[i] - 1]);
        }
    } catch (std::invalid_argument const& ex) {
        throw std::invalid_argument(fmt::format("Error: Column {} is not numeric\n", i + 1));
    }

    return BedRecord{chrom, start, end, data, col_count};
}

BedRecord parseBedRecord(const std::string& bedString, const int valInd1) {
    std::vector<std::string_view> res = split_bedstring(bedString);
    std::vector<std::string> fields(res.begin(), res.end());
    if (valInd1 >= fields.size()) {
        throw std::out_of_range(fmt::format("Error: Column {} does not exist\n", valInd1 + 1));
    }

    std::string chrom = fields[0];
    int start = std::stoi(fields[1]);
    int end = std::stoi(fields[2]);
    float val1;

    try {
        val1 = std::stof(fields[valInd1]);
    } catch (std::invalid_argument const& ex) {
        throw std::invalid_argument(fmt::format("Error: Column {} is not numeric\n", valInd1 + 1));
    }

    std::vector<float> data = {val1};
    return BedRecord{chrom, start, end, data, 1};
};

//' Parse a bed record into chr, start, and compressed beta and cov
//'
//' @param A line from a bed file
BedRecord parseBiscuitRecord(const std::string& bedString) {
    std::vector<std::string_view> res = split_bedstring(bedString);
    std::vector<std::string> fields(res.begin(), res.end());
    std::string chrom = fields[0];
    int start = std::stoi(fields[1]);
    int end = std::stoi(fields[2]);
    float beta = std::stof(fields[3]);
    float cov = std::stof(fields[4]);
    float m_count = std::round(cov * beta);

    std::vector<float> data = {beta, cov, m_count};
    return BedRecord{chrom, start, end, data, 3};
}

//' Parse a bed record into chr, start, and compressed beta and cov
//'
//' @param A line from a bed file
BedRecord parseCovRecord(const std::string& bedString) {
    std::vector<std::string_view> res = split_bedstring(bedString);
    std::vector<std::string> fields(res.begin(), res.end());

    if (fields.size() < 6) {
        Rcpp::stop("Found too few columns in data for bismark output - aligner may not be in bismark format");
    }

    std::string chrom = fields[0];
    int start = std::stoi(fields[1]);
    int end = std::stoi(fields[2]);
    float beta = std::stof(fields[3]) / 100;
    float m_count = std::stof(fields[4]);
    float u_count = std::stof(fields[5]);
    float cov = m_count + u_count;

    std::vector<float> data = {beta, cov, m_count};
    return BedRecord{chrom, start, end, data, 3};
}

/*   Archived   */
/*   make CpGID from chr and start   */
/*
std::string CpGID(BedLine& parsed_bedline) {

    std::stringstream cpg_stream;

    cpg_stream << parsed_bedline.chr << ":" << parsed_bedline.start;
    return cpg_stream.str();
}
*/
