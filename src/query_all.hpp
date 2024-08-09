#ifndef QUERY_ALL_H
#define QUERY_ALL_H

#if defined __cplusplus

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <cstdio>
#include <filesystem>
#include "query.hpp"
#include "parsers.hpp"
#include "decoders.hpp"
#include <unordered_map>
#include "../inst/include/khashl.h"
#include "../inst/include/iscream_types.h"

typedef struct CpG {
    int chr, start;
} CpG;

static kh_inline khint_t kh_eq_cpg(CpG a, CpG b) {
    return (a.chr == b.chr) & (a.start == b.start);
}

static kh_inline khint_t kh_hash_cpg(CpG cpg) {
    khint_t to_hash; 

    // Create input value for hashing function
    // Pull off last 16 bits of chr
    uint16_t chr = cpg.chr & 0xffff;
    // Pull off last 16 bits of start
    khint_t start = cpg.start & 0xffff;
    to_hash = (chr <<16) + start;

    return kh_hash_uint32(to_hash);
}

#ifndef __MAP_INIT
#define __MAP_INIT
KHASHL_MAP_INIT(static, khmap_t, khmap, CpG, int, kh_hash_cpg, kh_eq_cpg);
#endif /* ifndef __MAP_INIT */

template <class Mat>
class QueryAll {

private:

    std::unordered_map<std::string, int> chr_map;
    std::unordered_map<int, std::string> chr_rev_map;
    khmap_t *cpg_map;

    bool is_merged;
    int n_intervals, n_cpgs, chr_id, n_samples;
    Rcpp::CharacterVector seqnames, sample_names;
    Rcpp::IntegerVector start;

public:

    QueryAll();
    QueryAll(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions, const bool bismark, const bool merged, const bool sparse, const int nthreads);
    void populate_matrix(RegionQuery& query, int& col_n, const bool bismark);
    void print_mat(std::vector<std::vector<int>>& matrix, const std::string& matrix_name);
    void print_QueryAll();

    Mat cov_mat, m_mat;
    Rcpp::List assays;
    Rcpp::List wrap() {
        return Rcpp::List::create(
            Rcpp::_("M") = assays["M"],
            Rcpp::_("Cov") = assays["Cov"],
            Rcpp::_("pos") = start,
            Rcpp::_("chr") = seqnames,
            Rcpp::_("sampleNames") = sample_names
        );
    }
};

template <class Mat>
QueryAll<Mat>::QueryAll() {
    n_intervals = 0;
    n_cpgs = 0;
}

template <class Mat>
QueryAll<Mat>::QueryAll(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions, const bool bismark, const bool merged, const bool sparse, const int nthreads) {
    n_cpgs = 0;
    chr_id = 0;
    n_samples = bedfile_vec.size();
    n_intervals = regions.size();
    sample_names = Rcpp::wrap(bedfile_vec);
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

    for (int i = 0; i < sample_names.size(); i++) {
        std::filesystem::path sample_path = bedfile_vec[i];
        sample_names[i] = sample_path.extension() == ".gz" ? sample_path.stem().stem().string() : sample_path.stem().string();
    }

    if (sparse) {
        Rcpp::S4 cov_rmat = Rcpp::wrap(cov_mat);
        Rcpp::S4 M_rmat = Rcpp::wrap(m_mat);
        cov_rmat.slot("Dimnames") = Rcpp::List::create(rownames, sample_names);
        M_rmat.slot("Dimnames") = Rcpp::List::create(rownames, sample_names);
        assays = Rcpp::List::create(
            Rcpp::_["Cov"] = cov_rmat,
            Rcpp::_["M"] = M_rmat
        );
    } else {
        Rcpp::NumericMatrix cov_rmat = Rcpp::wrap(cov_mat);
        Rcpp::NumericMatrix M_rmat = Rcpp::wrap(m_mat);
        Rcpp::colnames(cov_rmat) = sample_names;
        Rcpp::colnames(M_rmat) = sample_names;
        Rcpp::rownames(cov_rmat) = rownames;
        Rcpp::rownames(M_rmat) = rownames;

        assays = Rcpp::List::create(
            Rcpp::_["Cov"] = cov_rmat,
            Rcpp::_["M"] = M_rmat
        );
    }

}

template <class Mat>
void QueryAll<Mat>::populate_matrix(RegionQuery& query, int& col_n, const bool bismark) {

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

template <class Mat>
void QueryAll<Mat>::print_mat(std::vector<std::vector<int>>& matrix, const std::string& matrix_name) {
    Rprintf("%s\n", matrix_name.c_str());
    for (int j = 0; j < matrix[0].size(); j++) {
        Rprintf("[%d, ] ", j);
        for (int i = 0; i < matrix.size(); i++) {
            Rprintf(" %d", matrix[i][j]);
        }
        Rprintf("\n");
    }
}

#endif /* __cplusplus */

#endif /* ifndef QUERY_ALL_H */
