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

//' Sum CpGs M values and coverage
void aggregate(RegionQuery& interval, float& total_beta, float& total_cov, bool mval) {

    for (std::string& each_cpg : interval.cpgs_in_interval) {
        BedLine parsed_cpg = parseBEDRecord(each_cpg);
        total_beta += mval ? (int) std::round(parsed_cpg.cov * parsed_cpg.beta) : parsed_cpg.beta;
        total_cov += parsed_cpg.cov;
    }
}

//' Aggregate CpGs within features
//' @param bedfiles A vector of bedfile paths
//' @param regions A vector of genomic regions
//' @export
// [[Rcpp::export]]
void agg_cpgs_file(std::vector<std::string>& bedfiles, std::vector<std::string>& regions) {

    printf("Aggregating %zu regions from %zu bedfiles\n", regions.size(), bedfiles.size());

    FILE *scmet_matrix;
    scmet_matrix = std::fopen("scmet2.tsv", "w");
    std::vector<RegionQuery> cpgs_in_file(0);

    for (std::string& bedfile_name : bedfiles) {

        std::filesystem::path bed_path = bedfile_name;
        cpgs_in_file = query_intervals(bedfile_name.c_str(), regions);

        for (RegionQuery interval : cpgs_in_file) {
            int total_m = 0;
            int total_cov = 0;
            aggregate(interval, total_m, total_cov);
            fprintf(scmet_matrix, "%s\t%s\t%d\t%d\n", interval.interval_str.c_str(), bed_path.stem().stem().c_str(), total_cov, total_m);
        }
    }
    fclose(scmet_matrix);
}

//' Aggregate CpGs within features
//' @param bedfiles A vector of bedfile paths
//' @param regions A vector of genomic regions
//' @param region_rownames Whether to set rownames to the regions
//' @param nthreads Number of cores to use. See details.
//'
//' @details
//' The optimal number of threads depends on the number of bedfiles.
//'
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame cpg_apply(std::vector<std::string>& bedfiles, Rcpp::CharacterVector& regions, bool calculate_m, bool region_rownames = false, int nthreads = 1) {
    printf("Aggregating %zu regions from %zu bedfiles\n", regions.size(), bedfiles.size());
    ssize_t rowsize = bedfiles.size() * regions.size();
    Rcpp::CharacterVector feature_col(rowsize);
    Rcpp::CharacterVector cell(rowsize);
    Rcpp::NumericVector total_reads(rowsize);
    Rcpp::NumericVector me_reads(rowsize);

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
        for (RegionQuery interval : cpgs_in_file) {
            float mut_m_val = 0;
            float mut_cov_val = 0;
            aggregate(interval, mut_m_val, mut_cov_val, calculate_m);

            feature_col[row_count] = interval.interval_str;
            cell[row_count] = bed_path.stem().stem().c_str();
            total_reads[row_count] = mut_total_cov;
            me_reads[row_count] = mut_total_m;
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
        Rcpp::Named(m_beta_colname) = me_reads
    );

    if (region_rownames) {
        result.attr("row.names") = feature_col;
    }

    return result;
}
