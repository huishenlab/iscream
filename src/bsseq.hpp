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
#include "../inst/include/iscream_types.h"

class BS {

private:

    typedef std::unordered_map<std::string, int> CpGMap;
    CpGMap cpg_map;

    int n_intervals;
    int n_cpgs;
    int n_samples;

public:

    BS();
    BS(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions);
    void populate_arma_cols(RegionQuery& query, int& row_n);
    void print_mat(std::vector<std::vector<int>>& matrix, const std::string& matrix_name);
    void print_BS();
    const int size() const;

    arma::umat cov_mat;
    arma::umat m_mat;
    Rcpp::List assays;
};

#endif /* __cplusplus */

#endif /* ifndef BSSEQ_H */

