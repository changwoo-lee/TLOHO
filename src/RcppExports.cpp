// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cholesky_diagonalcpp
arma::mat cholesky_diagonalcpp(arma::mat& R, arma::vec& diagadd);
RcppExport SEXP _TLOHO_cholesky_diagonalcpp(SEXP RSEXP, SEXP diagaddSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type diagadd(diagaddSEXP);
    rcpp_result_gen = Rcpp::wrap(cholesky_diagonalcpp(R, diagadd));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TLOHO_cholesky_diagonalcpp", (DL_FUNC) &_TLOHO_cholesky_diagonalcpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_TLOHO(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
