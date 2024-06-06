#include <Rcpp.h>

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

//' Get the numebr of available threads
//' @param verbose Whether to be verbose on available omp threads
//' @export
//' @examples
//' get_omp_threads(verbose = TRUE)
// [[Rcpp::export]]
int get_omp_threads(bool verbose) {
    if (verbose) {
      #ifndef _OPENMP
        Rprintf("This installation of data.table has not been compiled with OpenMP support.\n");
      #else
        Rprintf("  OpenMP version (_OPENMP)       %d\n", _OPENMP);
      #endif
      Rprintf("  omp_get_num_procs()            %d\n", omp_get_num_procs());
      Rprintf("  omp_get_thread_limit()         %d\n", omp_get_thread_limit());
      Rprintf("  omp_get_max_threads()          %d\n", omp_get_max_threads());
    }
    return omp_get_max_threads();
}
