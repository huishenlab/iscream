#ifndef QUERY_ALL_H
#define QUERY_ALL_H

#if defined __cplusplus

#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(stringfish)]]
#include <sf_external.h>

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <cstdio>
#include <filesystem>
#include "query.hpp"
#include "parsers.hpp"
#include "log.hpp"
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
    int n_intervals, n_cpgs, chr_id, n_samples, resize_count;
    Rcpp::CharacterVector sample_names;
    SEXP seqnames;
    Rcpp::IntegerVector start;

public:

    QueryAll();
    QueryAll(
        std::vector<std::string>& bedfile_vec,
        std::vector<std::string>& regions,
        const bool bismark,
        const bool merged,
        const bool sparse,
        const int prealloc,
        const int nthreads
    );
    void populate_matrix(RegionQuery& query, int& col_n, const bool bismark);
    void print_QueryAll();

    Mat bitmat;
    Rcpp::List assays;
    Rcpp::List wrap() {
        UNPROTECT(1); // seqnames
        return Rcpp::List::create(
            Rcpp::_("M") = assays["M"],
            Rcpp::_("packed") = assays["packed"],
            Rcpp::_("pos") = start,
            Rcpp::_("chr") = seqnames,
            Rcpp::_("sampleNames") = sample_names
        );
    }
};

#endif /* __cplusplus */

#endif /* ifndef QUERY_ALL_H */
