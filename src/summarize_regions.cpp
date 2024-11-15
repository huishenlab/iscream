// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <cstdio>
#include <cmath>
#include <filesystem>
#include "query.hpp"
#include "log.hpp"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
#include <spdlog/fmt/bundled/ranges.h>
// TODO: switch to ranges.h once spdlog 0.0.19 is available for ubuntu

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
    CPG_COUNT,
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
    {"cpg_count", CPG_COUNT},
    {"min", MIN},
    {"max", MAX},
    {"range", RANGE},
};

// colnames for the output DataFrame
typedef struct {
    std::string cov_name;
    std::string mval_name;
} Colnames;

// Container for computed summary value vectors
typedef struct {
    StatFunction fun;
    std::vector<double> coverage;
    std::vector<double> meth;
    Colnames colnames;
} ComputedVec;

// Input data vectors
typedef struct {
    arma::vec covs;
    arma::vec mvals;
} DataVec;


// Produce vectors of the input coverage and beta/M values
DataVec make_data_vec(const std::vector<std::string> cpgs, const bool mval, const bool bismark) {

    arma::vec covs(cpgs.size());
    arma::vec mvals(cpgs.size());

    for (int i = 0; i < cpgs.size(); i++) {
        BedLine parsed_cpg = bismark ? parseCovRecord(cpgs[i]) : parseBEDRecord(cpgs[i]);
        covs[i] = parsed_cpg.cov;
        mvals[i] = mval ? parsed_cpg.m_count : parsed_cpg.beta;
    }
    return DataVec{covs, mvals};
}

// Produce vectors pairs and colnames of computed summaries that will go into the result DataFrame
ComputedVec init_computed_vec(const int rowsize, const std::string fun, const bool mval) {
    std::string prefix = mval ? "M." : "beta.";
    return ComputedVec {
        str_to_enum[fun],
        std::vector<double>(rowsize, -99),
        std::vector<double>(rowsize, -99),
        Colnames {"coverage." + fun, prefix + fun}
    };
}

// Create the necessary vector pairs for requested summary functions
std::vector<ComputedVec> init_result_cols(const int rowsize, const std::vector<std::string> fun_vec, const bool mval) {
    std::vector<ComputedVec> vecs(fun_vec.size());
    for (int i = 0; i < fun_vec.size(); i++) {
        vecs[i] = init_computed_vec(rowsize, fun_vec[i], mval);
    }
    return vecs;
}

// Return {computed coverage, computed_mval}
std::tuple<double, double> compute_vecs(const StatFunction func, const DataVec& data_vec) {
    switch (func) {
        case MEAN:
            return {arma::mean(data_vec.covs), arma::mean(data_vec.mvals)};
        case MEDIAN:
            return {arma::median(data_vec.covs), arma::median(data_vec.mvals)};
        case STDDEV:
            return {arma::stddev(data_vec.covs), arma::stddev(data_vec.mvals)};
        case VARIANCE:
            return {arma::var(data_vec.covs), arma::var(data_vec.mvals)};
        case CPG_COUNT:
            return {data_vec.covs.size(), data_vec.mvals.size()};
        case MIN:
            return {arma::min(data_vec.covs), arma::min(data_vec.mvals)};
        case MAX:
            return {arma::max(data_vec.covs), arma::max(data_vec.mvals)};
        case RANGE:
            return {arma::range(data_vec.covs), arma::range(data_vec.mvals)};
        default: // using sum as default since R CMD check won't allow a switch without default
            return {arma::sum(data_vec.covs), arma::sum(data_vec.mvals)};
    }
}


//' Apply a function over CpGs within features
//'
//' This function should be called from `summarize_regions()` since there are few
//' sanity checks on the C++ side.
//' @param bedfiles A vector of bedfile paths
//' @param regions A vector of genomic regions
//' @param fun_vec Vector of the armadillo-supported stats functions to apply over the
//' CpGs in the ' regions: `"sum"`, `"mean"`, `"median"`, `"stddev"`,
//' `"variance"` "`cpg_count`", `"min"`,`"max"`, and `"range"`.
//' @param mval Calculates M values when TRUE, use beta values when FALSE
//' @param bismark If the input is in the bismark column format instead of BISCUIT
//' @param region_rownames Whether to set rownames to the regions strings. Not
//' necessary if your regions vector is unnamed. If its names, then the "Feature"
//' column is set to the names and the rownames are set to the regions string
//' @param nthreads Number of cores to use. See details.
//'
//' @details
//' The optimal number of threads depends on the number of bedfiles, but is set
//' to half the available OpenMP cores. See `?get_threads` for more details. It
//' can be manaully set with `set_threads()`.
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame Cpp_summarize_regions(
    const std::vector<std::string>& bedfiles,
    const Rcpp::CharacterVector& regions,
    const std::vector<std::string>& fun_vec,
    const bool mval,
    const bool bismark,
    const bool region_rownames = false,
    const int& nthreads = 1
) {

    // LOGGER
    setup_logger("iscream::summarize_regions");

    std::string fun_label = fmt::format("{}", fmt::join(fun_vec, ", "));
    spdlog::info("Summarizing {} regions from {} bedfiles", regions.size(), bedfiles.size());
    spdlog::info("using {}", fun_label.c_str());

    std::vector<std::string> regions_vec = Rcpp::as<std::vector<std::string>>(regions);
    ssize_t rowsize = bedfiles.size() * regions.size();

    spdlog::stopwatch sw;

    Rcpp::CharacterVector feature_col(rowsize, Rcpp::CharacterVector::get_na());
    Rcpp::CharacterVector sample(rowsize, Rcpp::CharacterVector::get_na());
    std::vector<ComputedVec> computed_vecs = init_result_cols(rowsize, fun_vec, mval);
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
                DataVec data_vec = make_data_vec(interval.cpgs_in_interval, mval, bismark);

                for (ComputedVec& vec : computed_vecs) {
                    std::tie(vec.coverage[row_count], vec.meth[row_count]) = compute_vecs(vec.fun, data_vec);
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
        Rcpp::Named("Feature") = (regions.hasAttribute("names") ? (Rcpp::CharacterVector) regions.names() : feature_col),
        Rcpp::Named("Sample") = sample
    );

    for (ComputedVec vec : computed_vecs) {
        result.push_back(vec.coverage, vec.colnames.cov_name);
        result.push_back(vec.meth, vec.colnames.mval_name);
    }
    spdlog::debug("Created DataFrame in {} s", sw);

    if (region_rownames) {
        result.attr("row.names") = feature_col;
        spdlog::debug("Rownames set");
    }

    return result;
}

