#include <Rcpp.h>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include "query.hpp"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    float computed_beta_me;
    float computed_cov;
} ComputedCpG;

// Sum CpGs M values and coverage
ComputedCpG aggregate(
    const RegionQuery& interval,
    const bool mval
) {

    float total_beta = 0;
    float total_cov = 0;
    for (const std::string& each_cpg : interval.cpgs_in_interval) {
        BedLine parsed_cpg = parseBEDRecord(each_cpg);
        total_beta += mval ? (int) std::round(parsed_cpg.cov * parsed_cpg.beta) : parsed_cpg.beta;
        total_cov += parsed_cpg.cov;
    }

    return ComputedCpG{total_beta, total_cov};
}

// Get mean of betas and coverage
ComputedCpG mean(
    const RegionQuery& interval,
    const bool mval
) {

    int n_cpg = interval.cpgs_in_interval.size();

    ComputedCpG cpg = aggregate(interval, mval);

    cpg.computed_beta_me = n_cpg == 0 ? 0.0 : cpg.computed_beta_me / n_cpg;
    cpg.computed_cov = n_cpg == 0 ? 0 : cpg.computed_cov / n_cpg;

    return cpg;
}

//' Apply a function over CpGs within features
//'
//' This function should be called from `region_map()` since there are few
//' sanity checks on the C++ side.
//' @param bedfiles A vector of bedfile paths
//' @param regions A vector of genomic regions
//' @param fun One of the supported functions to apply over the CpGs in the
//' regions: `"aggregate"`, `"average"`.
//' @param mval Calculates M values when TRUE, use beta values when FALSE
//' @param region_rownames Whether to set rownames to the regions strings
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
    const std::string& fun,
    const bool mval,
    const bool region_rownames = false,
    const int& nthreads = 1
) {

    Rprintf("Aggregating %zu regions from %zu bedfiles\n", regions.size(), bedfiles.size());
    ssize_t rowsize = bedfiles.size() * regions.size();
    Rcpp::CharacterVector feature_col(rowsize);
    Rcpp::CharacterVector cell(rowsize);
    Rcpp::NumericVector total_reads(rowsize);
    Rcpp::NumericVector beta_me_reads(rowsize);

    auto f = aggregate;
    if (fun == "average") {
        f = mean;
    }

    std::vector<std::string> regions_vec = Rcpp::as<std::vector<std::string>>(regions);
    Progress bar(bedfiles.size(), true);
    int completed_beds = 0;
#if defined(_OPENMP)
    #pragma omp parallel for num_threads(nthreads)
#endif
    for (int bedfile_n = 0; bedfile_n < bedfiles.size(); bedfile_n++) {
        std::vector<RegionQuery> cpgs_in_file(0);
        std::string bedfile_name = bedfiles[bedfile_n];
        std::filesystem::path bed_path = bedfile_name;
        cpgs_in_file = query_intervals(bedfile_name.c_str(), regions_vec);
        int empty_cpg_count = 0;

        int row_count = bedfile_n * regions_vec.size();
        std::string bedfile_prefix = bed_path.stem().stem();
        for (RegionQuery interval : cpgs_in_file) {
            ComputedCpG agg_cpg = f(interval, mval);

            feature_col[row_count] = interval.interval_str;
            cell[row_count] = bedfile_prefix.c_str();
            total_reads[row_count] = agg_cpg.computed_cov;
            beta_me_reads[row_count] = agg_cpg.computed_beta_me;

            row_count++;
            empty_cpg_count += interval.cpgs_in_interval.size();
        }
        // TODO: thread-safe way to warn when no cpgs are found in interval.
        // Lots of warnings from multiple threads cause stack overflow
        bar.increment();
    }

    Rcpp::String m_beta_colname = mval ? "M" : "beta";
    Rcpp::DataFrame result = Rcpp::DataFrame::create(
        Rcpp::Named("Feature") = (regions.hasAttribute("names") ? (Rcpp::CharacterVector) regions.names() : feature_col),
        Rcpp::Named("Cell") = cell,
        Rcpp::Named("coverage") = total_reads,
        Rcpp::Named(m_beta_colname) = beta_me_reads
    );

    if (region_rownames) {
        result.attr("row.names") = feature_col;
    }

    return result;
}
