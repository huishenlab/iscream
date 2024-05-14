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

BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type bedfile_vec(bedfile_vecSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type regions(regionsSEXP);
END_RCPP
}
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type bedfiles(bedfilesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type regions(regionsSEXP);
    Rcpp::traits::input_parameter< bool >::type region_rownames(region_rownamesSEXP);
    rcpp_result_gen = Rcpp::wrap(agg_cpgs_df(bedfiles, regions, region_rownames));
    return rcpp_result_gen;
END_RCPP
}
// decode_beta
int decode_beta(uint32_t encoded);
RcppExport SEXP _iscream_decode_beta(SEXP encodedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint32_t >::type encoded(encodedSEXP);
    rcpp_result_gen = Rcpp::wrap(decode_beta(encoded));
    return rcpp_result_gen;
END_RCPP
}
// decode_cov
int decode_cov(uint32_t encoded);
RcppExport SEXP _iscream_decode_cov(SEXP encodedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint32_t >::type encoded(encodedSEXP);
    rcpp_result_gen = Rcpp::wrap(decode_cov(encoded));
    return rcpp_result_gen;
END_RCPP
}
// decode_m
int decode_m(uint32_t encoded);
RcppExport SEXP _iscream_decode_m(SEXP encodedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< uint32_t >::type encoded(encodedSEXP);
    rcpp_result_gen = Rcpp::wrap(decode_m(encoded));
    return rcpp_result_gen;
END_RCPP
}
// vdecoder
IntegerVector vdecoder(IntegerVector& encoded, int measure);
RcppExport SEXP _iscream_vdecoder(SEXP encodedSEXP, SEXP measureSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type encoded(encodedSEXP);
    Rcpp::traits::input_parameter< int >::type measure(measureSEXP);
    rcpp_result_gen = Rcpp::wrap(vdecoder(encoded, measure));
    return rcpp_result_gen;
END_RCPP
}
// vdouble_decoder
DoubleVector vdouble_decoder(DoubleVector& encoded, int measure);
RcppExport SEXP _iscream_vdouble_decoder(SEXP encodedSEXP, SEXP measureSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DoubleVector& >::type encoded(encodedSEXP);
    Rcpp::traits::input_parameter< int >::type measure(measureSEXP);
    rcpp_result_gen = Rcpp::wrap(vdouble_decoder(encoded, measure));
    return rcpp_result_gen;
END_RCPP
}
// vencoder
IntegerVector vencoder(NumericVector beta_col, IntegerVector cov_col);
RcppExport SEXP _iscream_vencoder(SEXP beta_colSEXP, SEXP cov_colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_col(beta_colSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type cov_col(cov_colSEXP);
    rcpp_result_gen = Rcpp::wrap(vencoder(beta_col, cov_col));
    return rcpp_result_gen;
END_RCPP
}
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
// query_interval
std::vector<std::vector<std::string>> query_interval(std::vector<std::string>& bedfiles, std::string& region);
RcppExport SEXP _iscream_query_interval(SEXP bedfilesSEXP, SEXP regionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type bedfiles(bedfilesSEXP);
    Rcpp::traits::input_parameter< std::string& >::type region(regionSEXP);
    rcpp_result_gen = Rcpp::wrap(query_interval(bedfiles, region));
    return rcpp_result_gen;
END_RCPP
}
// cpg_apply
Rcpp::DataFrame cpg_apply(std::vector<std::string>& bedfiles, Rcpp::CharacterVector& regions, std::string fun, bool mval, bool region_rownames, int nthreads);
RcppExport SEXP _iscream_cpg_apply(SEXP bedfilesSEXP, SEXP regionsSEXP, SEXP funSEXP, SEXP mvalSEXP, SEXP region_rownamesSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type bedfiles(bedfilesSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector& >::type regions(regionsSEXP);
    Rcpp::traits::input_parameter< std::string >::type fun(funSEXP);
    Rcpp::traits::input_parameter< bool >::type mval(mvalSEXP);
    Rcpp::traits::input_parameter< bool >::type region_rownames(region_rownamesSEXP);
    Rcpp::traits::input_parameter< int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cpg_apply(bedfiles, regions, fun, mval, region_rownames, nthreads));
    return rcpp_result_gen;
END_RCPP
}
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_iscream_decode_beta", (DL_FUNC) &_iscream_decode_beta, 1},
    {"_iscream_decode_cov", (DL_FUNC) &_iscream_decode_cov, 1},
    {"_iscream_decode_m", (DL_FUNC) &_iscream_decode_m, 1},
    {"_iscream_vdecoder", (DL_FUNC) &_iscream_vdecoder, 2},
    {"_iscream_vdouble_decoder", (DL_FUNC) &_iscream_vdouble_decoder, 2},
    {"_iscream_vencoder", (DL_FUNC) &_iscream_vencoder, 2},
    {"_iscream_get_omp_threads", (DL_FUNC) &_iscream_get_omp_threads, 1},
    {"_iscream_query_interval", (DL_FUNC) &_iscream_query_interval, 2},
    {"_iscream_cpg_apply", (DL_FUNC) &_iscream_cpg_apply, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_iscream(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
