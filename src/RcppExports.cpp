// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/iscream_types.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// agg_cpgs_file
void agg_cpgs_file(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions);
RcppExport SEXP _iscream_agg_cpgs_file(SEXP bedfile_vecSEXP, SEXP regionsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type bedfile_vec(bedfile_vecSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type regions(regionsSEXP);
    agg_cpgs_file(bedfile_vec, regions);
    return R_NilValue;
END_RCPP
}
// agg_cpgs_df
Rcpp::DataFrame agg_cpgs_df(std::vector<std::string>& bedfile_vec, std::vector<std::string>& regions);
RcppExport SEXP _iscream_agg_cpgs_df(SEXP bedfile_vecSEXP, SEXP regionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type bedfile_vec(bedfile_vecSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type regions(regionsSEXP);
    rcpp_result_gen = Rcpp::wrap(agg_cpgs_df(bedfile_vec, regions));
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
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_iscream_agg_cpgs_file", (DL_FUNC) &_iscream_agg_cpgs_file, 2},
    {"_iscream_agg_cpgs_df", (DL_FUNC) &_iscream_agg_cpgs_df, 2},
    {"_iscream_decode_beta", (DL_FUNC) &_iscream_decode_beta, 1},
    {"_iscream_decode_cov", (DL_FUNC) &_iscream_decode_cov, 1},
    {"_iscream_decode_m", (DL_FUNC) &_iscream_decode_m, 1},
    {"_iscream_vdecoder", (DL_FUNC) &_iscream_vdecoder, 2},
    {"_iscream_vdouble_decoder", (DL_FUNC) &_iscream_vdouble_decoder, 2},
    {"_iscream_vencoder", (DL_FUNC) &_iscream_vencoder, 2},
    {"_iscream_query_interval", (DL_FUNC) &_iscream_query_interval, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_iscream(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
