#ifndef TREE_H
#define TREE_H

#if defined __cplusplus

#include "../inst/include/iscream_types.h"
#include <Rcpp.h>
#include "parsers.hpp"
#include "query.hpp"
#include <unordered_map>

typedef std::vector<std::vector<EncodedBedLine>> ReadMatrix;

class Tree {

private:

    struct CpG {
        CpG();
        CpG(int sample_number, int encoded_value);
        uint16_t sample;
        uint32_t encoded;
    };


    typedef std::unordered_map<std::string, std::vector<CpG>*> CpGMap;
    struct Interval {
        Interval(std::string& region, std::vector<std::string>& bedfile_vec);
        std::string interval_str;
        CpGMap cpg_map;
    };

    std::vector<Interval> intervals;
    int n_intervals;

public:
    Tree();
    Tree(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions);
    const int size() const;
    void printTree(const std::string& m_matrix = "M.mtx", const std::string cov_matrix = "cov.mtx");
};

#endif /* __cplusplus */

#endif /* ifndef TREE_H */
