// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <cstdio>
#include <cmath>
#include <filesystem>
#include "query.hpp"
#include "log.hpp"
#include "../inst/include/iscream_types.h"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <spdlog/fmt/bundled/ranges.h>

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

enum StatFunction {
    SUM,
    MEAN,
    MEDIAN,
    STDDEV,
    VARIANCE,
    COUNT,
    MIN,
    MAX,
    RANGE
};

// get Function enum from input 'fun' argument
std::unordered_map<std::string, StatFunction> str_to_enum {
    {"sum", SUM},
    {"mean", MEAN},
    {"median", MEDIAN},
    {"stddev", STDDEV},
    {"variance", VARIANCE},
    {"count", COUNT},
    {"min", MIN},
    {"max", MAX},
    {"range", RANGE},
};

// Input data vectors with length equal to the number of input functions
typedef std::vector<arma::vec> InputCols;

// Container for computed summary value vectors
typedef struct {
    std::string fun_string;
    StatFunction fun;
    std::vector<std::vector<double>> computed_vecs;
} ComputedFunVecs;

InputCols make_BSdata_vecs(
    const std::vector<std::string>& cpgs,
    const ssize_t rowsize,
    const std::vector<int> col_indices,
    const std::vector<std::string>& fun_vec,
    const BSType bstype,
    const bool mval
) {

    // only two arma vectors needed: [0] is cov and [1] is beta/M
    std::vector<arma::vec> input_data(2, arma::vec(cpgs.size()));
    BedRecord parsed_cpg;
    for (int i = 0; i < cpgs.size(); i++) {
        switch (bstype) {
            case BISMARK:
                parsed_cpg = parseCovRecord(cpgs[i]);
                break;
            default:
                parsed_cpg = parseBiscuitRecord(cpgs[i]);
        }
        input_data[0][i] = parsed_cpg.data[1];
        input_data[1][i] = mval ? parsed_cpg.data[2] : parsed_cpg.data[0];
    }
    return input_data;
}

InputCols make_data_vecs(
    const std::vector<std::string>& cpgs,
    const ssize_t rowsize,
    const std::vector<int> col_indices,
    const std::vector<std::string>& fun_vec,
    const bool mval
) {

    std::vector<arma::vec> input_data(col_indices.size(), arma::vec(cpgs.size()));

    for (int i = 0; i < cpgs.size(); i++) {
        BedRecord parsed_cpg = parseBedRecord(cpgs[i], col_indices);
        for (int dp = 0; dp < parsed_cpg.size; dp++) {
            input_data[dp][i] = parsed_cpg.data[dp];
        }
    }
    return input_data;
}

ComputedFunVecs init_computed_vec(
    const ssize_t rowsize,
    const ssize_t colsize,
    const std::string& fun,
    const bool mval
) {
    std::vector<std::vector<double>> computed_vecs(colsize, std::vector<double>(rowsize, -99));
    return ComputedFunVecs {
        fun,
        str_to_enum[fun],
        computed_vecs
    };
}

std::vector<ComputedFunVecs> init_result_cols(
    const ssize_t rowsize,
    const ssize_t colsize,
    const std::vector<std::string>& fun_vec,
    const bool mval
) {

    std::vector<ComputedFunVecs> vecs(fun_vec.size());
    for (int fun_ind = 0; fun_ind < fun_vec.size(); fun_ind++) {
        vecs[fun_ind] = init_computed_vec(rowsize, colsize, fun_vec[fun_ind], mval);
    }
    return vecs;
}

// Return {computed coverage, computed_mval}
double summarize(const StatFunction func, const arma::vec& data_vec) {
    switch (func) {
        case MEAN:
            return arma::mean(data_vec);
        case MEDIAN:
            return arma::median(data_vec);
        case STDDEV:
            return arma::stddev(data_vec);
        case VARIANCE:
            return arma::var(data_vec);
        case COUNT:
            return data_vec.size();
        case MIN:
            return arma::min(data_vec);
        case MAX:
            return arma::max(data_vec);
        case RANGE:
            return arma::range(data_vec);
        default: // using sum as default since R CMD check won't allow a switch without default
            return arma::sum(data_vec);
    }
}


//' Apply a function over BED file records within genomic features
//'
//' This function should be called from `summarize_regions()` since there are few
//' sanity checks on the C++ side.
//' @param bedfiles A vector of bedfile paths
//' @param regions A vector of genomic regions
//' @param fun_vec Vector of the armadillo-supported stats functions to apply over the
//' CpGs in the ' regions: `"sum"`, `"mean"`, `"median"`, `"stddev"`,
//' `"variance"` "`count`", `"min"`,`"max"`, and `"range"`.
//' @param col_indices A vector of genomic regions
//' @param col_names A vector of genomic regions
//' @param mval Calculates M values when TRUE, use beta values when FALSE
//' @param region_rownames Whether to set rownames to the regions strings. Not
//' necessary if your regions vector is unnamed. If its names, then the "feature"
//' column is set to the names and the rownames are set to the regions string
//' @param nthreads Number of cores to use. See details.
//'
//' @details
//' The optimal number of threads depends on the number of bedfiles, but is set
//' to half the available OpenMP cores. See `?get_threads` for more details. It
//' can be manaully set with `set_threads()`.
//'
//' @returns A summary data.frame
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame Cpp_summarize_regions(
    const std::vector<std::string>& bedfiles,
    const Rcpp::CharacterVector& regions,
    const std::vector<std::string>& fun_vec,
    const std::vector<int>& col_indices,
    const std::vector<std::string>& col_names,
    const std::string& aligner,
    const bool mval = false,
    const bool region_rownames = false,
    const int nthreads = 1
) {

    // LOGGER
    setup_logger("iscream::summarize_regions");

    std::string fun_label = fmt::format("{}", fmt::join(fun_vec, ", "));
    spdlog::info("Summarizing {} regions from {} bedfiles", regions.size(), bedfiles.size());
    spdlog::info("using {}", fun_label.c_str());
    spdlog::info(
        "with columns {} as {}",
        fmt::format("{}", fmt::join(col_indices, ", ")),
        fmt::format("{}", fmt::join(col_names, ", "))
    );

    BSType bstype = get_BSType(aligner);
    std::vector<std::string> regions_vec = Rcpp::as<std::vector<std::string>>(regions);
    ssize_t rowsize = bedfiles.size() * regions.size();

    spdlog::stopwatch sw;

    Rcpp::CharacterVector feature_col(rowsize, Rcpp::CharacterVector::get_na());
    Rcpp::CharacterVector sample(rowsize, Rcpp::CharacterVector::get_na());
    std::vector<ComputedFunVecs> computed_vecs = init_result_cols(rowsize, col_indices.size(), fun_vec, mval);
    spdlog::debug("Created vectors for DataFrame with {} rows in {} s", rowsize, sw);

    sw.reset();

    Progress bar(bedfiles.size(), true);
    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int bedfile_n = 0; bedfile_n < bedfiles.size(); bedfile_n++) {
        if ( !Progress::check_abort() ) {
            std::vector<RegionQuery> cpgs_in_file(0);
            std::string bedfile_name = bedfiles[bedfile_n];
            std::filesystem::path bed_path = bedfile_name;
            cpgs_in_file = query_intervals(bedfile_name.c_str(), regions_vec);
            int empty_cpg_count = 0;

            int row_count = bedfile_n * regions_vec.size();
            std::string bedfile_prefix = bed_path.stem().stem();
            spdlog::debug("Got {} as sample name from {}", bedfile_prefix, bedfile_name);

            for (RegionQuery interval : cpgs_in_file) {
                spdlog::debug("Got {} CpGs from {}", interval.cpgs_in_interval.size(), bedfile_name);
                feature_col[row_count] = interval.interval_str;
                sample[row_count] = bedfile_prefix.c_str();

                if (interval.cpgs_in_interval.size() == 0) {
                    row_count++;
                    continue;
                }

                InputCols input_data;
                switch(bstype) {
                    case BISCUIT: case BISMARK:
                        input_data = make_BSdata_vecs(interval.cpgs_in_interval, rowsize, col_indices, fun_vec, bstype, mval);
                        break;
                    default:
                        input_data = make_data_vecs(interval.cpgs_in_interval, rowsize, col_indices, fun_vec, mval);
                }

                for (ComputedFunVecs& vecs : computed_vecs) {
                    for (int i = 0; i < input_data.size(); i++) {
                        vecs.computed_vecs[i][row_count] = summarize(vecs.fun, input_data[i]);
                    }
                }

                row_count++;
                empty_cpg_count += interval.cpgs_in_interval.size();
            }
            // TODO: thread-safe way to warn when no cpgs are found in interval.
            // Lots of warnings from multiple threads cause stack overflow
            if (spdlog::get_level() == spdlog::level::info) bar.increment();
        }
    }

    spdlog::debug("Populated columns in {} s", sw);
    sw.reset();

    Rcpp::DataFrame result = Rcpp::DataFrame::create(
        Rcpp::Named("feature") = (regions.hasAttribute("names") ? (Rcpp::CharacterVector) regions.names() : feature_col),
        Rcpp::Named("file") = sample
    );

    for (int i = 0; i < computed_vecs.size(); i++) {
        ComputedFunVecs vecs = computed_vecs[i];
        for (int j = 0; j < vecs.computed_vecs.size(); j++) {
            std::vector<double> vec = vecs.computed_vecs[j];
            result.push_back(vec, col_names[j] + "." + vecs.fun_string);
        }
    }
    spdlog::debug("Created DataFrame in {} s", sw);

    if (region_rownames) {
        result.attr("row.names") = feature_col;
        spdlog::debug("Rownames set");
    }

    return result;
}

