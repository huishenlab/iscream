// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/iscream_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_omp_threads
int get_omp_threads(bool verbose);
RcppExport SEXP _iscream_get_omp_threads(SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(get_omp_threads(verbose));
    return rcpp_result_gen;
END_RCPP
}
// setup_logger
void setup_logger(std::string logname);
RcppExport SEXP _iscream_setup_logger(SEXP lognameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type logname(lognameSEXP);
    setup_logger(logname);
    return R_NilValue;
END_RCPP
}
// Cpp_set_log_level
void Cpp_set_log_level(const std::string& name);
RcppExport SEXP _iscream_Cpp_set_log_level(SEXP nameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type name(nameSEXP);
    Cpp_set_log_level(name);
    return R_NilValue;
END_RCPP
}
// get_log_level
std::string get_log_level();
RcppExport SEXP _iscream_get_log_level() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(get_log_level());
    return rcpp_result_gen;
END_RCPP
}
// Cpp_query_chroms
std::set<std::string> Cpp_query_chroms(const std::vector<std::string>& bedfile_vec, const int nthreads);
RcppExport SEXP _iscream_Cpp_query_chroms(SEXP bedfile_vecSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type bedfile_vec(bedfile_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Cpp_query_chroms(bedfile_vec, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// Cpp_query_interval
Rcpp::CharacterVector Cpp_query_interval(const std::string& bedfile, const std::vector<std::string>& regions);
RcppExport SEXP _iscream_Cpp_query_interval(SEXP bedfileSEXP, SEXP regionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type bedfile(bedfileSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type regions(regionsSEXP);
    rcpp_result_gen = Rcpp::wrap(Cpp_query_interval(bedfile, regions));
    return rcpp_result_gen;
END_RCPP
}
// scan_tabix
Rcpp::List scan_tabix(const std::string& bedfile, const std::vector<std::string>& regions);
RcppExport SEXP _iscream_scan_tabix(SEXP bedfileSEXP, SEXP regionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type bedfile(bedfileSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type regions(regionsSEXP);
    rcpp_result_gen = Rcpp::wrap(scan_tabix(bedfile, regions));
    return rcpp_result_gen;
END_RCPP
}
// Cpp_query_all
Rcpp::List Cpp_query_all(std::vector<std::string>& bedfiles, std::vector<std::string>& regions, const bool bismark, const bool merged, const bool sparse, const int prealloc, const int nthreads);
RcppExport SEXP _iscream_Cpp_query_all(SEXP bedfilesSEXP, SEXP regionsSEXP, SEXP bismarkSEXP, SEXP mergedSEXP, SEXP sparseSEXP, SEXP preallocSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type bedfiles(bedfilesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type regions(regionsSEXP);
    Rcpp::traits::input_parameter< const bool >::type bismark(bismarkSEXP);
    Rcpp::traits::input_parameter< const bool >::type merged(mergedSEXP);
    Rcpp::traits::input_parameter< const bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< const int >::type prealloc(preallocSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Cpp_query_all(bedfiles, regions, bismark, merged, sparse, prealloc, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// Cpp_summarize_regions
Rcpp::DataFrame Cpp_summarize_regions(const std::vector<std::string>& bedfiles, const Rcpp::CharacterVector& regions, const std::vector<std::string>& fun_vec, const bool mval, const bool bismark, const bool region_rownames, const int& nthreads);
RcppExport SEXP _iscream_Cpp_summarize_regions(SEXP bedfilesSEXP, SEXP regionsSEXP, SEXP fun_vecSEXP, SEXP mvalSEXP, SEXP bismarkSEXP, SEXP region_rownamesSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type bedfiles(bedfilesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type regions(regionsSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type fun_vec(fun_vecSEXP);
    Rcpp::traits::input_parameter< const bool >::type mval(mvalSEXP);
    Rcpp::traits::input_parameter< const bool >::type bismark(bismarkSEXP);
    Rcpp::traits::input_parameter< const bool >::type region_rownames(region_rownamesSEXP);
    Rcpp::traits::input_parameter< const int& >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(Cpp_summarize_regions(bedfiles, regions, fun_vec, mval, bismark, region_rownames, nthreads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_iscream_get_omp_threads", (DL_FUNC) &_iscream_get_omp_threads, 1},
    {"_iscream_setup_logger", (DL_FUNC) &_iscream_setup_logger, 1},
    {"_iscream_Cpp_set_log_level", (DL_FUNC) &_iscream_Cpp_set_log_level, 1},
    {"_iscream_get_log_level", (DL_FUNC) &_iscream_get_log_level, 0},
    {"_iscream_Cpp_query_chroms", (DL_FUNC) &_iscream_Cpp_query_chroms, 2},
    {"_iscream_Cpp_query_interval", (DL_FUNC) &_iscream_Cpp_query_interval, 2},
    {"_iscream_scan_tabix", (DL_FUNC) &_iscream_scan_tabix, 2},
    {"_iscream_Cpp_query_all", (DL_FUNC) &_iscream_Cpp_query_all, 7},
    {"_iscream_Cpp_summarize_regions", (DL_FUNC) &_iscream_Cpp_summarize_regions, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_iscream(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
