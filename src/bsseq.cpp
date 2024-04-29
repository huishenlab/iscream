#include <cstdio>
#include "bsseq.hpp"
#include <filesystem>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

BS::BS() {
    n_intervals = 0;
    n_cpgs = 0;
}

BS::BS(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions) {
    n_cpgs = 0;
    n_samples = bedfile_vec.size();
    n_intervals = regions.size();
    sample_names = bedfile_vec;

    cov_mat.resize(5, bedfile_vec.size());
    m_mat.resize(5, bedfile_vec.size());

    printf("Querying %zu regions from %zu bedfiles\n", regions.size(), bedfile_vec.size());
    Progress bar(bedfile_vec.size(), true);

    for (int bedfile_n = 0; bedfile_n < bedfile_vec.size(); bedfile_n++) {
        MultiRegionQuery cpgs_in_file = query_intervals(bedfile_vec[bedfile_n].c_str(), regions);
        for (RegionQuery cpgs_in_interval : cpgs_in_file) {
            populate_arma_cols(cpgs_in_interval, bedfile_n);
        }
        bar.increment();
    }

    if (cov_mat.n_rows > cpg_map.size()) {
        printf("Correcting matrix size\n");
        int diff_rows = cov_mat.n_rows - cpg_map.size();
        cov_mat.resize(cpg_map.size(), bedfile_vec.size());
        m_mat.resize(cpg_map.size(), bedfile_vec.size());
    }

    Rcpp::NumericMatrix cov_rmat = Rcpp::wrap(cov_mat);
    Rcpp::NumericMatrix M_rmat = Rcpp::wrap(m_mat);
    for (int i = 0; i < sample_names.size(); i++) {
        std::filesystem::path sample_path = sample_names[i];
        sample_names[i] = sample_path.extension() == ".gz" ? sample_path.stem().stem().string() : sample_path.stem().string();
    }

    Rcpp::colnames(cov_rmat) = Rcpp::wrap(sample_names);
    Rcpp::colnames(M_rmat) = Rcpp::wrap(sample_names);
    assays = Rcpp::List::create(
        Rcpp::_["Cov"] = cov_rmat,
        Rcpp::_["M"] = M_rmat
    );
}

void BS::populate_arma_cols(RegionQuery& query, int& col_n) {

    for (std::string cpg_string : query.cpgs_in_interval) {
        BedLine parsed_bedline = parseBEDRecord(cpg_string);
        std::string cpg_id = CpGID(parsed_bedline);
        if (!cpg_map.count(cpg_id)) {
            n_cpgs++;
            cpg_map.insert({cpg_id, n_cpgs});
            chrs.push_back(parsed_bedline.chr);
            starts.push_back(parsed_bedline.start);
        }

        if (cov_mat.n_rows < cpg_map.size()) {
            int extra_rows = cpg_map.size() - cov_mat.n_rows;
            cov_mat.resize(cov_mat.n_rows + extra_rows, cov_mat.n_cols);
            m_mat.resize(m_mat.n_rows + extra_rows, cov_mat.n_cols);
        }

        cov_mat(cpg_map[cpg_id] - 1, col_n) = parsed_bedline.cov;
        m_mat(cpg_map[cpg_id] - 1, col_n) =  (int) std::round(parsed_bedline.cov * parsed_bedline.beta);
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
    printf("Cov\n");
    cov_mat.print();
    printf("M\n");
    m_mat.print();
}

const int BS::size() const {
    return n_intervals;
}
