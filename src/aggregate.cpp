#include <Rcpp.h>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include "query.hpp"

void aggregate(RegionQuery& interval, int& total_m, int& total_cov) {

    for (std::string& each_cpg : interval.cpgs_in_interval) {
        BedLine parsed_cpg = parseBEDRecord(each_cpg);
        total_m += (int) std::round(parsed_cpg.cov * parsed_cpg.beta);
        total_cov += parsed_cpg.cov;
    }
}

//' Aggregate CpGs within features
//' @param bedfiles A vector of bedfile paths
//' @param regions A vector of genomic regions
//' @export
// [[Rcpp::export]]
void agg_cpgs_file(std::vector<std::string>& bedfiles, std::vector<std::string>& regions) {
    printf("n_bedfiles: %zu\n", bedfiles.size());
    printf("n_intervals: %zu\n", regions.size());

    FILE *scmet_matrix;
    scmet_matrix = std::fopen("scmet2.tsv", "w");
    std::vector<RegionQuery> cpgs_in_file(0);

    for (std::string& bedfile_name : bedfiles) {

        std::filesystem::path bed_path = bedfile_name;
        printf("Querying %s\n", bedfile_name.c_str());
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
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame agg_cpgs_df(std::vector<std::string>& bedfiles, std::vector<std::string>& regions) {
    printf("n_bedfiles: %zu\n", bedfiles.size());
    printf("n_intervals: %zu\n", regions.size());

    ssize_t rowsize = bedfiles.size() * regions.size();
    std::vector<std::string> feature_col(rowsize);
    std::vector<std::string> cell(rowsize);
    std::vector<int> total_reads(rowsize);
    std::vector<int> me_reads(rowsize);

    std::vector<RegionQuery> cpgs_in_file(0);
    int row_count = 0;

    for (std::string& bedfile_name : bedfiles) {

        std::filesystem::path bed_path = bedfile_name;
        printf("Querying %s\n", bedfile_name.c_str());
        cpgs_in_file = query_intervals(bedfile_name.c_str(), regions);

        for (RegionQuery interval : cpgs_in_file) {
            int mut_total_m = 0;
            int mut_total_cov = 0;
            aggregate(interval, mut_total_m, mut_total_cov);
            feature_col[row_count] = interval.interval_str;
            cell[row_count] = bed_path.stem().stem();
            total_reads[row_count] = mut_total_cov;
            me_reads[row_count] = mut_total_m;
            row_count++;
        }
    }

    Rcpp::DataFrame result = Rcpp::DataFrame::create(
        Rcpp::Named("Feature") = feature_col,
        Rcpp::Named("Cell") = cell,
        Rcpp::Named("total_reads") = total_reads,
        Rcpp::Named("met_reads") = me_reads
    );
    return result;
}
