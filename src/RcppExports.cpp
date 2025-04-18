// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// allele_diff_strings
std::vector<std::string> allele_diff_strings(std::vector<std::string> germs, int X, Rcpp::Nullable<Rcpp::CharacterVector> non_mismatch_chars_nullable);
RcppExport SEXP _piglet_allele_diff_strings(SEXP germsSEXP, SEXP XSEXP, SEXP non_mismatch_chars_nullableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type germs(germsSEXP);
    Rcpp::traits::input_parameter< int >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type non_mismatch_chars_nullable(non_mismatch_chars_nullableSEXP);
    rcpp_result_gen = Rcpp::wrap(allele_diff_strings(germs, X, non_mismatch_chars_nullable));
    return rcpp_result_gen;
END_RCPP
}
// allele_diff_indices
std::vector<int> allele_diff_indices(std::vector<std::string> germs, int X, Rcpp::Nullable<Rcpp::CharacterVector> non_mismatch_chars_nullable);
RcppExport SEXP _piglet_allele_diff_indices(SEXP germsSEXP, SEXP XSEXP, SEXP non_mismatch_chars_nullableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type germs(germsSEXP);
    Rcpp::traits::input_parameter< int >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type non_mismatch_chars_nullable(non_mismatch_chars_nullableSEXP);
    rcpp_result_gen = Rcpp::wrap(allele_diff_indices(germs, X, non_mismatch_chars_nullable));
    return rcpp_result_gen;
END_RCPP
}
// allele_diff_indices_parallel
Rcpp::RObject allele_diff_indices_parallel(std::vector<std::string> germs, std::vector<std::string> inputs, int X, bool parallel, bool return_count);
RcppExport SEXP _piglet_allele_diff_indices_parallel(SEXP germsSEXP, SEXP inputsSEXP, SEXP XSEXP, SEXP parallelSEXP, SEXP return_countSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type germs(germsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type inputs(inputsSEXP);
    Rcpp::traits::input_parameter< int >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< bool >::type return_count(return_countSEXP);
    rcpp_result_gen = Rcpp::wrap(allele_diff_indices_parallel(germs, inputs, X, parallel, return_count));
    return rcpp_result_gen;
END_RCPP
}
// allele_diff_indices_parallel2
Rcpp::RObject allele_diff_indices_parallel2(std::vector<std::string> germs, std::vector<std::string> inputs, int X, bool parallel, bool return_count, Rcpp::Nullable<Rcpp::CharacterVector> non_mismatch_chars_nullable);
RcppExport SEXP _piglet_allele_diff_indices_parallel2(SEXP germsSEXP, SEXP inputsSEXP, SEXP XSEXP, SEXP parallelSEXP, SEXP return_countSEXP, SEXP non_mismatch_chars_nullableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type germs(germsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type inputs(inputsSEXP);
    Rcpp::traits::input_parameter< int >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< bool >::type return_count(return_countSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::CharacterVector> >::type non_mismatch_chars_nullable(non_mismatch_chars_nullableSEXP);
    rcpp_result_gen = Rcpp::wrap(allele_diff_indices_parallel2(germs, inputs, X, parallel, return_count, non_mismatch_chars_nullable));
    return rcpp_result_gen;
END_RCPP
}
// insert_gaps2_vec
std::vector<std::string> insert_gaps2_vec(const std::vector<std::string>& gapped, const std::vector<std::string>& ungapped, bool parallel);
RcppExport SEXP _piglet_insert_gaps2_vec(SEXP gappedSEXP, SEXP ungappedSEXP, SEXP parallelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type gapped(gappedSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type ungapped(ungappedSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    rcpp_result_gen = Rcpp::wrap(insert_gaps2_vec(gapped, ungapped, parallel));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_piglet_allele_diff_strings", (DL_FUNC) &_piglet_allele_diff_strings, 3},
    {"_piglet_allele_diff_indices", (DL_FUNC) &_piglet_allele_diff_indices, 3},
    {"_piglet_allele_diff_indices_parallel", (DL_FUNC) &_piglet_allele_diff_indices_parallel, 5},
    {"_piglet_allele_diff_indices_parallel2", (DL_FUNC) &_piglet_allele_diff_indices_parallel2, 6},
    {"_piglet_insert_gaps2_vec", (DL_FUNC) &_piglet_insert_gaps2_vec, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_piglet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
