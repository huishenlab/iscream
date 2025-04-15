#ifndef ISCREAM_TYPES_H
#define ISCREAM_TYPES_H
#if defined __cplusplus

#include <string>

// define aligner types
enum BSType {
    BISMARK,
    BISCUIT,
    GENERAL
};

// get alinger
inline BSType get_BSType(const std::string& aligner) {
    if (aligner == "biscuit") {
        return BISCUIT;
    } else if (aligner == "bismark" || aligner == "bsbolt") {
        return BISMARK;
    } else {
        return GENERAL;
    }
}

#endif /* __cplusplus */
#endif /* ifndef QUERY_ALL_H */
