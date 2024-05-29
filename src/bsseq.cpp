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
    chr_id = 0;
    n_samples = bedfile_vec.size();
    n_intervals = regions.size();
    sample_names = bedfile_vec;
    cpg_map = khmap_init();

    cov_mat.resize(5, bedfile_vec.size());
    m_mat.resize(5, bedfile_vec.size());

    printf("Querying %zu regions from %zu bedfiles\n", regions.size(), bedfile_vec.size());
    Progress bar(bedfile_vec.size(), true); 

    for (int bedfile_n = 0; bedfile_n < bedfile_vec.size(); bedfile_n++) {
        MultiRegionQuery cpgs_in_file = query_intervals(bedfile_vec[bedfile_n].c_str(), regions);
        for (RegionQuery cpgs_in_interval : cpgs_in_file) {
            populate_matrix(cpgs_in_interval, bedfile_n);
        }
         bar.increment(); 
    }

    int mapsize = kh_size(cpg_map);
    if (cov_mat.n_rows > mapsize) {
        printf("Correcting matrix size\n");
        int diff_rows = cov_mat.n_rows - mapsize;
        cov_mat.resize(mapsize, bedfile_vec.size());
        m_mat.resize(mapsize, bedfile_vec.size());
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

void BS::populate_matrix(RegionQuery& query, int& col_n) {

    std::vector<BedLine> lines;
    std::vector<CpG> ids;
    khmap_m_resize(cpg_map, query.cpgs_in_interval.size());
    for (std::string cpg_string : query.cpgs_in_interval) {

        BedLine parsed_bedline = parseBEDRecord(cpg_string);
        lines.push_back(parsed_bedline);

        unsigned int i = 0;
        if (!chr_map.count(parsed_bedline.chr)) chr_map.insert({parsed_bedline.chr, ++chr_id});

        CpG cpg = CpG{chr_map[parsed_bedline.chr], parsed_bedline.start};

        ids.push_back(cpg);

        khint_t insert_b;
        int absent;
        insert_b = khmap_put(cpg_map, cpg, &absent);
        if (absent) {
            n_cpgs++;
            kh_val(cpg_map, insert_b) = n_cpgs;
            chrs.push_back(parsed_bedline.chr);
            starts.push_back(parsed_bedline.start);
        }
    }

    int mapsize = kh_size(cpg_map);
    if (cov_mat.n_rows < mapsize) {
        int extra_rows = (mapsize - cov_mat.n_rows) * 1000;
        cov_mat.resize(cov_mat.n_rows + extra_rows, cov_mat.n_cols);
        m_mat.resize(m_mat.n_rows + extra_rows, cov_mat.n_cols);
    }

    khint_t retrieve_b;
    int idx;
    for (size_t i = 0; i < lines.size(); i++) {
        retrieve_b = khmap_get(cpg_map, ids[i]);
        idx = kh_val(cpg_map, retrieve_b);
        cov_mat(idx - 1, col_n) = lines[i].cov;
        m_mat(idx - 1, col_n) =  (int) std::round(lines[i].cov * lines[i].beta);
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
