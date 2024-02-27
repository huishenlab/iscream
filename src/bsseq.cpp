#include <cstdio>
#include "bsseq.hpp"
#include "../inst/include/indicators.hpp"

BS::BS() {
    n_intervals = 0;
    n_cpgs = 0;
}

BS::BS(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions) {
    n_cpgs = 0;
    n_samples = bedfile_vec.size();
    n_intervals = regions.size();
    sample_names = bedfile_vec;

    indicators::ProgressBar bar {
        indicators::option::BarWidth{50},
        indicators::option::Start{"["},
        indicators::option::Fill{"°"},
        indicators::option::Lead{" "},
        indicators::option::Remainder{" "},
        indicators::option::End{"]"},
        indicators::option::ShowPercentage{true},
        indicators::option::ForegroundColor{indicators::Color::cyan},
        indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}
    };

    cov_mat.resize(5, bedfile_vec.size());
    m_mat.resize(5, bedfile_vec.size());

    float completed_beds = 0;
    printf("Querying %zu regions from %zu bedfiles\n", regions.size(), bedfile_vec.size());
    bar.set_progress(0);
    for (int bedfile_n = 0; bedfile_n < bedfile_vec.size(); bedfile_n++) {
        MultiRegionQuery cpgs_in_file = query_intervals(bedfile_vec[bedfile_n].c_str(), regions);
        for (RegionQuery cpgs_in_interval : cpgs_in_file) {
            populate_arma_cols(cpgs_in_interval, bedfile_n);
        }
        completed_beds++;
        if (completed_beds < bedfile_vec.size()) {
            bar.set_option(indicators::option::PostfixText{bedfile_vec[bedfile_n]});
            bar.set_progress((int) std::round(completed_beds / bedfile_vec.size() * 100));
        } else {
            bar.set_option(indicators::option::PostfixText{"Done: " + std::to_string(cpg_map.size()) + " CpGs found!"});
            bar.set_progress(100);
        }
    }

    if (cov_mat.n_rows > cpg_map.size()) {
        printf("Correcting matrix size\n");
        int diff_rows = cov_mat.n_rows - cpg_map.size();
        cov_mat.resize(cpg_map.size(), bedfile_vec.size());
        m_mat.resize(cpg_map.size(), bedfile_vec.size());
    }

    Rcpp::NumericMatrix cov_rmat = Rcpp::wrap(cov_mat);
    Rcpp::NumericMatrix M_rmat = Rcpp::wrap(m_mat);
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
