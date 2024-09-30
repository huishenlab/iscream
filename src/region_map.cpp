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

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

enum Function {
    SUM,
    MEAN,
    MEDIAN,
    STDDEV,
    VARIANCE,
    RANGE
};

// get Function enum from input 'fun' argument
std::unordered_map<std::string, Function> str_to_enum {
    {"sum", SUM},
    {"mean", MEAN},
    {"median", MEDIAN},
    {"stddev", STDDEV},
    {"variance", VARIANCE},
    {"range", RANGE},
};

// colnames for the output DataFrame
typedef struct {
    std::string cov_name;
    std::string mval_name;
} Colnames;

// Container for computed summary value vectors
typedef struct {
    Function func;
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
DataVec makevec(const std::vector<std::string> cpgs, const bool mval, const bool bismark) {

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
ComputedVec get_vec(const int rowsize, const std::string fun, const bool mval) {
    std::string prefix = mval ? "M." : "beta.";
    return ComputedVec {
        str_to_enum[fun],
        std::vector<double>(rowsize, -99),
        std::vector<double>(rowsize, -99),
        Colnames {"coverage." + fun, prefix + fun}
    };
}

// Create the necessary vector pairs for requested summary functions
std::vector<ComputedVec> get_vectors(const int rowsize, const std::vector<std::string> funcs, const bool mval) {
    std::vector<ComputedVec> vecs(funcs.size());
    for (int i = 0; i < funcs.size(); i++) {
        vecs[i] = get_vec(rowsize, funcs[i], mval);
    }
    return vecs;
}

//' Apply a function over CpGs within features
//'
//' This function should be called from `region_map()` since there are few
//' sanity checks on the C++ side.
//' @param bedfiles A vector of bedfile paths
//' @param regions A vector of genomic regions
//' @param fun One of the armadillo-supported stats functions to apply over the
//' CpGs in the ' regions: `"sum"`, `"mean"`, `"median"`, `"stddev"`,
//' `"variance"`, `"range"`.
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
Rcpp::DataFrame Cpp_region_map(
    const std::vector<std::string>& bedfiles,
    const Rcpp::CharacterVector& regions,
    const std::vector<std::string>& funcs,
    const bool mval,
    const bool bismark,
    const bool region_rownames = false,
    const int& nthreads = 1
) {

    // LOGGER
    setup_logger("iscream::region_map");

    std::string func_label;
    for (const std::string& fun : funcs) func_label += fun + " ";

    spdlog::info("{} {} regions from {} bedfiles\n", func_label.c_str(), regions.size(), bedfiles.size());
    ssize_t rowsize = bedfiles.size() * regions.size();

    Rcpp::CharacterVector feature_col(rowsize, Rcpp::CharacterVector::get_na());
    Rcpp::CharacterVector cell(rowsize, Rcpp::CharacterVector::get_na());
    std::vector<ComputedVec> computed_vecs = get_vectors(rowsize, funcs, mval);
    spdlog::debug("Created vectors for DataFrame with {} rows in {} s", rowsize, sw);

    std::vector<std::string> regions_vec = Rcpp::as<std::vector<std::string>>(regions);
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
                cell[row_count] = bedfile_prefix.c_str();

                if (interval.cpgs_in_interval.size() == 0) {
                    row_count++;
                    continue;
                }
                DataVec data_vec = makevec(interval.cpgs_in_interval, mval, bismark);

                for (ComputedVec& vec : computed_vecs) {
                    switch (vec.func) {
                        case SUM:
                            vec.coverage[row_count] = arma::sum(data_vec.covs);
                            vec.meth[row_count] = arma::sum(data_vec.mvals);
                            break;
                        case MEAN:
                            vec.coverage[row_count] = arma::mean(data_vec.covs);
                            vec.meth[row_count] = arma::mean(data_vec.mvals);
                            break;
                        case MEDIAN:
                            vec.coverage[row_count] = arma::median(data_vec.covs);
                            vec.meth[row_count] = arma::median(data_vec.mvals);
                            break;
                        case STDDEV:
                            vec.coverage[row_count] = arma::stddev(data_vec.covs);
                            vec.meth[row_count] = arma::stddev(data_vec.mvals);
                            break;
                        case VARIANCE:
                            vec.coverage[row_count] = arma::var(data_vec.covs);
                            vec.meth[row_count] = arma::var(data_vec.mvals);
                            break;
                        case RANGE:
                            vec.coverage[row_count] = arma::range(data_vec.covs);
                            vec.meth[row_count] = arma::range(data_vec.mvals);
                            break;
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

    std::string m_beta_colname = mval ? "M" : "beta";
    spdlog::debug("Using {} as the methylation value name", m_beta_colname);
    Rcpp::DataFrame result = Rcpp::DataFrame::create(
        Rcpp::Named("Feature") = (regions.hasAttribute("names") ? (Rcpp::CharacterVector) regions.names() : feature_col),
        Rcpp::Named("Cell") = cell
    );

    for (ComputedVec vec : computed_vecs) {
        result.push_back(vec.coverage, vec.colnames.cov_name);
        result.push_back(vec.meth, vec.colnames.mval_name);
    }

    if (region_rownames) {
        result.attr("row.names") = feature_col;
        spdlog::debug("Rownames set");
    }

    return result;
}

