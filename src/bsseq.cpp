#include <cstdio>
#include "bsseq.hpp"
#include <filesystem>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

BS::BS() {
    n_intervals = 0;
    n_cpgs = 0;
}

BS::BS(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions, const bool bismark, const bool merged, const int nthreads) {
    n_cpgs = 0;
    chr_id = 0;
    n_samples = bedfile_vec.size();
    n_intervals = regions.size();
    sample_names = bedfile_vec;
    cpg_map = khmap_init();
    is_merged = merged;

    cov_mat.resize(5, bedfile_vec.size());
    m_mat.resize(5, bedfile_vec.size());

    Rprintf("Querying %zu regions from %zu bedfiles\n", regions.size(), bedfile_vec.size());
    Progress bar(bedfile_vec.size(), true); 

#if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthreads)
#endif
    for (int bedfile_n = 0; bedfile_n < bedfile_vec.size(); bedfile_n++) {
        MultiRegionQuery cpgs_in_file = query_intervals(bedfile_vec[bedfile_n].c_str(), regions);
        for (RegionQuery cpgs_in_interval : cpgs_in_file) {
            populate_matrix(cpgs_in_interval, bedfile_n, bismark);
        }
        bar.increment();
    }
    bar.cleanup();

    Rcpp::CharacterVector c(kh_size(cpg_map));
    Rcpp::IntegerVector s(kh_size(cpg_map));
    Rcpp::CharacterVector rownames(kh_size(cpg_map));

    khint_t iter;
    for (iter = 0; iter < kh_end(cpg_map); ++iter) {
        if (kh_exist(cpg_map, iter)) {
            CpG cpg = kh_key(cpg_map, iter);
            int row_idx = kh_val(cpg_map, iter) - 1;
            s[row_idx] = cpg.start;
            c[row_idx] = chr_rev_map[cpg.chr];
            std::stringstream cpgid_stream;
            cpgid_stream << chr_rev_map[cpg.chr] << ":" << cpg.start + 1;
            rownames[row_idx] = cpgid_stream.str();
        }
    }
    seqnames = c;
    start = s;

    int mapsize = kh_size(cpg_map);
    if (cov_mat.n_rows > mapsize) {
        Rprintf("Correcting matrix size\n");
        int diff_rows = cov_mat.n_rows - mapsize;
        cov_mat.resize(mapsize, bedfile_vec.size());
        m_mat.resize(mapsize, bedfile_vec.size());
    }

    Rcpp::NumericMatrix cov_rmat = Rcpp::wrap(cov_mat);
    Rcpp::NumericMatrix M_rmat = Rcpp::wrap(m_mat);
    Rcpp::rownames(cov_rmat) = rownames;
    Rcpp::rownames(M_rmat) = rownames;

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

void BS::populate_matrix(RegionQuery& query, int& col_n, const bool bismark) {

    std::vector<BedLine> lines;
    std::vector<CpG> ids;
    #pragma omp critical
    {
        khmap_m_resize(cpg_map, query.cpgs_in_interval.size());
    }
    for (std::string cpg_string : query.cpgs_in_interval) {

        BedLine parsed_bedline = bismark ? parseCovRecord(cpg_string) : parseBEDRecord(cpg_string);
        lines.push_back(parsed_bedline);

        #pragma omp critical
        {
            if (!chr_map.count(parsed_bedline.chr)) {
                chr_map.insert({parsed_bedline.chr, ++chr_id});
                chr_rev_map.insert({chr_id, parsed_bedline.chr});
            }
        }

        CpG cpg = CpG{chr_map[parsed_bedline.chr], parsed_bedline.start};

        ids.push_back(cpg);

        #pragma omp critical
        {
            khint_t insert_b;
            int absent;
            insert_b = khmap_put(cpg_map, cpg, &absent);
            if (absent) {
                n_cpgs++;
                kh_val(cpg_map, insert_b) = n_cpgs;
            }
        }
    }

    #pragma omp critical
    {
        int mapsize = kh_size(cpg_map);
        if (cov_mat.n_rows < mapsize) {
            int extra_rows = (mapsize - cov_mat.n_rows) * 1000;
            cov_mat.resize(cov_mat.n_rows + extra_rows, cov_mat.n_cols);
            m_mat.resize(m_mat.n_rows + extra_rows, cov_mat.n_cols);
        }
    }

    khint_t retrieve_b;
    int idx;
    for (size_t i = 0; i < lines.size(); i++) {
        retrieve_b = khmap_get(cpg_map, ids[i]);
        idx = kh_val(cpg_map, retrieve_b);
        cov_mat(idx - 1, col_n) = lines[i].cov;
        m_mat(idx - 1, col_n) =  lines[i].m_count;
    }
}

void BS::print_mat(std::vector<std::vector<int>>& matrix, const std::string& matrix_name) {
    Rprintf("%s\n", matrix_name.c_str());
    for (int j = 0; j < matrix[0].size(); j++) {
        Rprintf("[%d, ] ", j);
        for (int i = 0; i < matrix.size(); i++) {
            Rprintf(" %d", matrix[i][j]);
        }
        Rprintf("\n");
    }
}

void BS::print_BS() {
    Rprintf("Cov\n");
    cov_mat.print();
    Rprintf("M\n");
    m_mat.print();
}
