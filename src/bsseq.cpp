#include <cstdio>
#include "bsseq.hpp"

BS::BS() {
    n_intervals = 0;
    n_cpgs = 0;
}

BS::BS(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions) {

    n_intervals  = 0;
    n_cpgs = 0;
    n_samples = bedfile_vec.size();
    cov_matrix_cols.resize(bedfile_vec.size(), std::vector<int>(10));
    m_matrix_cols.resize(bedfile_vec.size(), std::vector<int>(10));

    printf("n_bedfiles: %zu\n", bedfile_vec.size());
    printf("n_intervals: %zu\n", regions.size());
    for (int bedfile_n = 0; bedfile_n < bedfile_vec.size(); bedfile_n++) {
        printf("file: %s\n", bedfile_vec[bedfile_n].c_str());
        MultiRegionQuery cpgs_in_file = query_intervals(bedfile_vec[bedfile_n].c_str(), regions);
        for (RegionQuery cpgs_in_interval : cpgs_in_file) {
            populate_matrix_cols(cpgs_in_interval, bedfile_n);
        }
    }
    if (cov_matrix_cols[0].size() > cpg_map.size()) {
        resize_mats(cpg_map.size());
    }
}

void BS::populate_matrix_cols(RegionQuery& query, int& col_n) {
    for (std::string cpg_string : query.cpgs_in_interval) {

        BedLine parsed_bedline = parseBEDRecord(cpg_string);
        std::string cpg_id = CpGID(parsed_bedline);
        if (!cpg_map.count(cpg_id)) {
            n_cpgs++;
            cpg_map.insert({cpg_id, n_cpgs});
        }

        if (cov_matrix_cols[col_n].size() < cpg_map.size()) {
            resize_mats(cov_matrix_cols[col_n].size() + query.cpgs_in_interval.size());
        }

        cov_matrix_cols[col_n][cpg_map[cpg_id] - 1] = parsed_bedline.cov;
        m_matrix_cols[col_n][cpg_map[cpg_id] - 1] = (int) std::round(parsed_bedline.cov * parsed_bedline.beta);
    }
}

void BS::resize_mats(const size_t newsize) {
    for (int col_n = 0; col_n < cov_matrix_cols.size(); col_n++) {
        cov_matrix_cols[col_n].resize(newsize);
        m_matrix_cols[col_n].resize(newsize);
    }
}


void BS::print_mat(std::vector<std::vector<int>>& matrix, const std::string& matrix_name) {
    printf("%s\n", matrix_name.c_str());
    for (int j = 0; j < matrix[0].size(); j++) {
        printf("[%d, ] ", j);
        for (int i = 0; i < matrix.size(); i++) {
            printf(" %d", matrix[i][j]);
        }
        printf("\n");
    }
}

void BS::print_BS() {
    print_mat(cov_matrix_cols, "Cov");
    print_mat(m_matrix_cols, "M");
}
}

const int BS::size() const {
    return n_intervals;
}
