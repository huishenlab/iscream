#ifndef BSSEQ_H
#define BSSEQ_H

#if defined __cplusplus

#include <cmath>
#include <Rcpp.h>
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
    void populate_matrix_cols(RegionQuery& query, int& row_n);
    void resize_mats(const size_t newsize);
    void print_mat(std::vector<std::vector<int>>& matrix, const std::string& matrix_name);
    void print_BS();
    const int size() const;

    std::vector<std::vector<int>> cov_matrix_cols;
    std::vector<std::vector<int>> m_matrix_cols;
};

#endif /* __cplusplus */

#endif /* ifndef BSSEQ_H */

