#include <Rcpp.h>

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// Get the number of available threads from OpenMP.
//
// This queries the number of available threads usign OpenMP, but will not
// reliably provide an accurate available thread count. To get a more reliable
// count that accounts for environment variables and HPC schedulers, use
// get_threads()`. This function was pulled from
// github.com/rdatatable/data.table
// @param verbose Whether to be verbose on available omp threads
//
// [[Rcpp::export]]
int get_omp_threads(bool verbose) {
    if (verbose) {
      #ifndef _OPENMP
        Rprintf("This installation of iscream has not been compiled with OpenMP support.\n");
      #else
        Rprintf("  OpenMP version (_OPENMP)       %d\n", _OPENMP);
        Rprintf("  omp_get_num_procs()            %d\n", omp_get_num_procs());
        Rprintf("  omp_get_thread_limit()         %d\n", omp_get_thread_limit());
        Rprintf("  omp_get_max_threads()          %d\n", omp_get_max_threads());
        return omp_get_max_threads();
      #endif
    }
    return 1;
}
