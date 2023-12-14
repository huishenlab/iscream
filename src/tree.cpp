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

Tree::CpG::CpG() {
    sample = 0;
    encoded = 0;
}

Tree::CpG::CpG(int sample_number, int encoded_value) {
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
        for (int cpg_n = 0; cpg_n < cpgs_1region_1file.size(); cpg_n++) {

            BedLine parsed_bedline = parseBEDRecord(cpgs_1region_1file[cpg_n]);
            std::string cpg_id = CpGID(parsed_bedline);

            if (!cpg_map.count(cpg_id)) {
                std::vector<CpG>* cpg_data_vec = new std::vector<CpG>;
                cpg_map.insert({cpg_id, cpg_data_vec});
            }

            EncodedBedLine encoded_cpg_line = encodeBedRecord(cpgs_1region_1file[cpg_n]);
            CpG tmp_cpg = CpG(bedfile_n, encoded_cpg_line.encoded);

            cpg_map[cpg_id]->push_back(tmp_cpg);
        }
    }
}

const int Tree::size() const {
    return n_intervals;
}

void Tree::printTree(const std::string& m_matrix, const std::string cov_matrix) {

    printf("%s\n", "printing");
    FILE *m_matrix_io;
    FILE *cov_matrix_io;

    m_matrix_io = std::fopen(m_matrix.c_str(), "w");
    cov_matrix_io = std::fopen(cov_matrix.c_str(), "w");

    for (Interval i : intervals) {
        for (auto item : i.cpg_map) {
            for (CpG cpg : *item.second) {
                fprintf(m_matrix_io, "%s\t%d\t%d\n", item.first.c_str(), cpg.sample + 1, decode_m(cpg.encoded));
                fprintf(cov_matrix_io, "%s\t%d\t%d\n", item.first.c_str(), cpg.sample + 1, decode_cov(cpg.encoded));
            }
            delete(item.second);
        }
    }

    fclose(m_matrix_io);
    fclose(cov_matrix_io);
}
