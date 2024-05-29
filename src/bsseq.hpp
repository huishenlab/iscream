#ifndef BSSEQ_H
#define BSSEQ_H

#if defined __cplusplus

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include "query.hpp"
#include "parsers.hpp"
#include "decoders.hpp"
#include <unordered_map>
#include "../inst/include/khashl.h"
#include "../inst/include/iscream_types.h"

#ifndef __MAP_INIT
#define __MAP_INIT
KHASHL_MAP_INIT(static, khmap_t, khmap, char*, int, kh_hash_str, kh_eq_str);
#endif /* ifndef __MAP_INIT */

class BS {

private:

    khmap_t *cpg_map;

    int n_intervals;
    int n_cpgs;
    int n_samples;
    Rcpp::CharacterVector chrs;
    Rcpp::IntegerVector starts;
    std::vector<std::string> sample_names;

public:

    BS();
    BS(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions);
    void populate_matrix(RegionQuery& query, int& col_n);
    void print_mat(std::vector<std::vector<int>>& matrix, const std::string& matrix_name);
    void print_BS();

    arma::umat cov_mat;
    arma::umat m_mat;
    Rcpp::List assays;
    Rcpp::S4 wrap() {
        Rcpp::Function BSseq("BSseq", Rcpp::Environment::namespace_env("bsseq"));
        Rcpp::Function GRanges("GRanges", Rcpp::Environment::namespace_env("GenomicRanges"));
        Rcpp::Function IRanges("IRanges", Rcpp::Environment::namespace_env("IRanges"));
        Rcpp::IntegerVector start = starts;
        Rcpp::IntegerVector end = starts;
        Rcpp::CharacterVector seqnames = chrs;
        return BSseq(
            Rcpp::_("M") = assays["M"],
            Rcpp::_("Cov") = assays["Cov"],
            Rcpp::_("gr") = GRanges(seqnames, IRanges(start, end)),
            Rcpp::_("sampleNames") = sample_names
        );
    }
};

#endif /* __cplusplus */

#endif /* ifndef BSSEQ_H */

