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
    struct DataPoint {
        DataPoint();
        DataPoint(int sample_number, int encoded_value);
        uint16_t sample;
        uint32_t encoded;
    };

    struct Interval {
        Interval(std::string& region, std::vector<std::string>& bedfile_vec);
        std::string interval_str;
        std::string chr;
        int start;
        int end;
        int n_cpgs;

        std::unordered_map<std::string, std::vector<DataPoint>*> cpg_map;
    };

    std::vector<Interval> intervals;
    int n_intervals;

public:
    Tree();
    Tree(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions);
    const int size() const;
    void printTree();
};

#endif /* __cplusplus */

#endif /* ifndef TREE_H */
