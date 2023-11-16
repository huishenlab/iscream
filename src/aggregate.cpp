#include <Rcpp.h>
#include <cstdio>
#include <cmath>
#include <filesystem>
#include "query.hpp"

//' @export
// [[Rcpp::export]]
void agg_cpgs_file(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions) {
    printf("n_bedfiles: %zu\n", bedfile_vec.size());
    printf("n_intervals: %zu\n", regions.size());

    FILE *scmet_matrix;
    scmet_matrix = std::fopen("scmet2.tsv", "w");
    std::vector<MultiRegionQuery> cpgs_in_file(0);

    for (std::string& bedfile_name : bedfile_vec) {

        std::filesystem::path bed_path = bedfile_name;
        cpgs_in_file = query_intervals(bedfile_name.c_str(), regions);

        for (MultiRegionQuery interval : cpgs_in_file) {
            int total_m = 0;
            int total_cov = 0;
            printf("%s\n", bedfile_name.c_str());
            for (std::string& each_cpg : interval.cpgs_in_interval) {
                BedLine parsed_cpg = parseBEDRecord(each_cpg);
                total_m += (int) std::round(parsed_cpg.cov * parsed_cpg.beta);
                total_cov += parsed_cpg.cov;
            }
            fprintf(scmet_matrix, "%s\t%s\t%d\t%d\n", interval.interval_str.c_str(), bed_path.stem().stem().c_str(), total_cov, total_m);
        }
    }
    fclose(scmet_matrix);
}

