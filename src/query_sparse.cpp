// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

#include "query.hpp"
#include "parsers.hpp"
#include "log.hpp"
#include <unordered_map>
#include "../inst/include/iscream_types.h"

//' @export
// [[Rcpp::export]]
void get_loci(
    const std::vector<std::string>& bedfiles,
    const Rcpp::CharacterVector& regions,
    const int nthreads = 1
) {
    std::vector<std::string> regions_vec = Rcpp::as<std::vector<std::string>>(regions);

    int count = 0;
    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(nthreads)
    #endif
    for (int bedfile_n = 0; bedfile_n < bedfiles.size(); bedfile_n++) {
        std::vector<RegionQuery> cpgs_in_file(0);
        std::string bedfile_name = bedfiles[bedfile_n];
        cpgs_in_file = query_intervals(bedfile_name.c_str(), regions_vec);

        for (RegionQuery interval : cpgs_in_file) {
            #pragma omp atomic
            count += interval.cpgs_in_interval.size();
        }
    }
    printf("%d\n", count);
}


