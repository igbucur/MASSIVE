// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// neg_log_prior_MR
double neg_log_prior_MR(List param_list, List prior_sd);
RcppExport SEXP _MASSIVE_neg_log_prior_MR(SEXP param_listSEXP, SEXP prior_sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type param_list(param_listSEXP);
    Rcpp::traits::input_parameter< List >::type prior_sd(prior_sdSEXP);
    rcpp_result_gen = Rcpp::wrap(neg_log_prior_MR(param_list, prior_sd));
    return rcpp_result_gen;
END_RCPP
}
// scaled_nl_posterior_MR
double scaled_nl_posterior_MR(unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, List param_list, List prior_sd, unsigned n);
RcppExport SEXP _MASSIVE_scaled_nl_posterior_MR(SEXP JSEXP, SEXP NSEXP, SEXP SSSEXP, SEXP sigma_GSEXP, SEXP param_listSEXP, SEXP prior_sdSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_G(sigma_GSEXP);
    Rcpp::traits::input_parameter< List >::type param_list(param_listSEXP);
    Rcpp::traits::input_parameter< List >::type prior_sd(prior_sdSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(scaled_nl_posterior_MR(J, N, SS, sigma_G, param_list, prior_sd, n));
    return rcpp_result_gen;
END_RCPP
}
// scaled_nl_gradient_MR
arma::vec scaled_nl_gradient_MR(unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, List param_list, List prior_sd, unsigned n);
RcppExport SEXP _MASSIVE_scaled_nl_gradient_MR(SEXP JSEXP, SEXP NSEXP, SEXP SSSEXP, SEXP sigma_GSEXP, SEXP param_listSEXP, SEXP prior_sdSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_G(sigma_GSEXP);
    Rcpp::traits::input_parameter< List >::type param_list(param_listSEXP);
    Rcpp::traits::input_parameter< List >::type prior_sd(prior_sdSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(scaled_nl_gradient_MR(J, N, SS, sigma_G, param_list, prior_sd, n));
    return rcpp_result_gen;
END_RCPP
}
// scaled_nl_hessian_MR
arma::mat scaled_nl_hessian_MR(unsigned J, unsigned N, arma::mat SS, arma::vec sigma_G, List param_list, List prior_sd, unsigned n);
RcppExport SEXP _MASSIVE_scaled_nl_hessian_MR(SEXP JSEXP, SEXP NSEXP, SEXP SSSEXP, SEXP sigma_GSEXP, SEXP param_listSEXP, SEXP prior_sdSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned >::type J(JSEXP);
    Rcpp::traits::input_parameter< unsigned >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type SS(SSSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_G(sigma_GSEXP);
    Rcpp::traits::input_parameter< List >::type param_list(param_listSEXP);
    Rcpp::traits::input_parameter< List >::type prior_sd(prior_sdSEXP);
    Rcpp::traits::input_parameter< unsigned >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(scaled_nl_hessian_MR(J, N, SS, sigma_G, param_list, prior_sd, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MASSIVE_neg_log_prior_MR", (DL_FUNC) &_MASSIVE_neg_log_prior_MR, 2},
    {"_MASSIVE_scaled_nl_posterior_MR", (DL_FUNC) &_MASSIVE_scaled_nl_posterior_MR, 7},
    {"_MASSIVE_scaled_nl_gradient_MR", (DL_FUNC) &_MASSIVE_scaled_nl_gradient_MR, 7},
    {"_MASSIVE_scaled_nl_hessian_MR", (DL_FUNC) &_MASSIVE_scaled_nl_hessian_MR, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_MASSIVE(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
