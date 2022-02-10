// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Sbruteforce
double Sbruteforce(NumericVector f, IntegerVector alpha);
RcppExport SEXP _numberofalleles_Sbruteforce(SEXP fSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type f(fSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(Sbruteforce(f, alpha));
    return rcpp_result_gen;
END_RCPP
}
// pr_sum
NumericVector pr_sum(IntegerVector x1, NumericVector fx1, IntegerVector x2, NumericVector fx2);
RcppExport SEXP _numberofalleles_pr_sum(SEXP x1SEXP, SEXP fx1SEXP, SEXP x2SEXP, SEXP fx2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fx1(fx1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fx2(fx2SEXP);
    rcpp_result_gen = Rcpp::wrap(pr_sum(x1, fx1, x2, fx2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_numberofalleles_Sbruteforce", (DL_FUNC) &_numberofalleles_Sbruteforce, 2},
    {"_numberofalleles_pr_sum", (DL_FUNC) &_numberofalleles_pr_sum, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_numberofalleles(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}