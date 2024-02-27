#include <Rcpp.h>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include "query.hpp"
#include "../inst/include/indicators.hpp"

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
//' @export
// [[Rcpp::export]]
Rcpp::DataFrame agg_cpgs_df(std::vector<std::string>& bedfiles, std::vector<std::string>& regions) {
    printf("Aggregating %zu regions from %zu bedfiles\n", regions.size(), bedfiles.size());

    ssize_t rowsize = bedfiles.size() * regions.size();
    Rcpp::CharacterVector feature_col(rowsize);
    Rcpp::CharacterVector cell(rowsize);
    Rcpp::NumericVector total_reads(rowsize);
    Rcpp::NumericVector me_reads(rowsize);

    std::vector<RegionQuery> cpgs_in_file(0);
    int row_count = 0;

    indicators::ProgressBar bar {
        indicators::option::BarWidth{50},
        indicators::option::Start{"["},
        indicators::option::Fill{"°"},
        indicators::option::Lead{" "},
        indicators::option::Remainder{" "},
        indicators::option::End{"]"},
        indicators::option::ShowPercentage{true},
        indicators::option::ForegroundColor{indicators::Color::cyan},
        indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}
    };

    float completed_beds = 0;
    for (std::string& bedfile_name : bedfiles) {

        std::filesystem::path bed_path = bedfile_name;
        std::vector<std::string> regions_vec = Rcpp::as<std::vector<std::string>>(regions);
        cpgs_in_file = query_intervals(bedfile_name.c_str(), regions_vec);

        for (RegionQuery interval : cpgs_in_file) {
            int mut_total_m = 0;
            int mut_total_cov = 0;
            aggregate(interval, mut_total_m, mut_total_cov);
            feature_col[row_count] = interval.interval_str;
            cell[row_count] = bed_path.stem().stem().c_str();
            total_reads[row_count] = mut_total_cov;
            me_reads[row_count] = mut_total_m;
            row_count++;
        }
        completed_beds++;
        if (completed_beds < bedfiles.size()) {
            bar.set_option(indicators::option::PostfixText{bedfile_name});
            bar.set_progress((int) std::round(completed_beds / bedfiles.size() * 100));
        } else {
            bar.set_option(indicators::option::PostfixText{"Done!"});
            bar.set_progress(100);
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
