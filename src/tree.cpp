#include <cstdio>
#include "tree.hpp"
#include "decoders.hpp"

Tree::Tree() {
    n_intervals = 0;
}

Tree::Tree(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions) {

    n_intervals  = 0;

    printf("n_bedfiles: %zu\n", bedfile_vec.size());
    printf("n_intervals: %zu\n", regions.size());
    for (std::string& region : regions) {
        printf("interval: %s\n", region.c_str());
        intervals.push_back(Interval(region, bedfile_vec));
        n_intervals++;
    }

}

Tree::DataPoint::DataPoint() {
    sample = 0;
    encoded = 0;
}


Tree::DataPoint::DataPoint(int sample_number, int encoded_value) {
    sample = sample_number;
    encoded = encoded_value;
}

Tree::Interval::Interval(std::string& line, std::vector<std::string>& bedfile_vec) {
    interval_str = line;
    std::istringstream ss(line);

    // get all cpgs in the 'line' region from every file
    for (int bedfile_n = 0; bedfile_n < bedfile_vec.size(); bedfile_n++) {
        printf("%d: %s\n", bedfile_n, bedfile_vec[bedfile_n].c_str());
        // get all cpgs from the 'line' region from bedfile_vec[i]
        std::vector<std::string> cpgs_1region_1file = query_interval(bedfile_vec[bedfile_n], line);

        // for each cpg in the line region from this file
        // add them to the hashmap with the encoded value and sample number
        for (int cpg = 0; cpg < cpgs_1region_1file.size(); cpg++) {

            std::string cpg_id = CpGID(cpgs_1region_1file[cpg]);

            if (!cpg_map.count(cpg_id)) {
                std::vector<DataPoint>* cpg_data_vec = new std::vector<DataPoint>;
                cpg_map.insert({cpg_id, cpg_data_vec});
            }

            EncodedBedLine encoded_cpg_line = encodeBedRecord(cpgs_1region_1file[cpg]);
            DataPoint tmp_dp = DataPoint(bedfile_n, encoded_cpg_line.encoded);

            cpg_map[cpg_id]->push_back(tmp_dp);
        }
    }
}

const int Tree::size() const {
    return n_intervals;
}

void Tree::printTree() {

    printf("%s\n", "printing");
    FILE *m_matrix;
    FILE *cov_matrix;

    m_matrix = std::fopen("M.mtx", "w");
    cov_matrix = std::fopen("cov.mtx", "w");

    for (Interval i : intervals) {
        for (auto item : i.cpg_map) {
            for (DataPoint dp : *item.second) {
                fprintf(m_matrix, "%s\t%d\t%d\n", item.first.c_str(), dp.sample + 1, decode_m(dp.encoded));
                fprintf(cov_matrix, "%s\t%d\t%d\n", item.first.c_str(), dp.sample + 1, decode_cov(dp.encoded));
            }
            delete(item.second);
        }
    }

    fclose(m_matrix);
    fclose(cov_matrix);
}
