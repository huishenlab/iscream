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
    Rcpp::CharacterVector chrs;
    Rcpp::IntegerVector starts;
    std::vector<std::string> sample_names;

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
            Rcpp::_("gr") = GRanges(seqnames, IRanges(start, end))
        );
    }

    Rcpp::S4 nobs_wrap() {

        // @assays
        Rcpp::S4 bss4(std::string("BSseq"));
        Rcpp::S4 simple_assays(std::string("SimpleAssays"));
        Rcpp::S4 data(std::string("SimpleList"));
        data.slot("listData") = assays;
        simple_assays.slot("data") = data;
        bss4.slot("assays") = simple_assays;

        // @colData
        Rcpp::S4 colData(std::string("DFrame"));
        colData.attr("nrows") = n_samples;
        colData.attr("rownames") = sample_names;
        bss4.slot("colData") = colData;

        // @elementMetadata
        Rcpp::S4 elementMetadata(std::string("DFrame"));
        elementMetadata.attr("nrows") = cpg_map.size();
        bss4.slot("elementMetadata") = elementMetadata;

        return bss4;
    }
};

#endif /* __cplusplus */

#endif /* ifndef BSSEQ_H */

