// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_facility_idx
List compute_facility_idx(arma::colvec& facility, int num_facility);
RcppExport SEXP _tdrecur_compute_facility_idx(SEXP facilitySEXP, SEXP num_facilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_facility_idx(facility, num_facility));
    return rcpp_result_gen;
END_RCPP
}
// update_beta
void update_beta(double l, double lchange, arma::colvec& L1, arma::mat& L2, arma::rowvec& Sm, arma::colvec& t_start, arma::colvec& t_end, arma::mat events_per_day_facility, arma::mat& z, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, arma::colvec& beta, arma::colvec& update, arma::colvec& update_rel, List facility_idx, arma::mat& hosp_begin, arma::colvec& max_d, bool parallel, const unsigned int& nthreads);
RcppExport SEXP _tdrecur_update_beta(SEXP lSEXP, SEXP lchangeSEXP, SEXP L1SEXP, SEXP L2SEXP, SEXP SmSEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP events_per_day_facilitySEXP, SEXP zSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP betaSEXP, SEXP updateSEXP, SEXP update_relSEXP, SEXP facility_idxSEXP, SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP parallelSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< double >::type lchange(lchangeSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type L2(L2SEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type Sm(SmSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type events_per_day_facility(events_per_day_facilitySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type update(updateSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type update_rel(update_relSEXP);
    Rcpp::traits::input_parameter< List >::type facility_idx(facility_idxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nthreads(nthreadsSEXP);
    update_beta(l, lchange, L1, L2, Sm, t_start, t_end, events_per_day_facility, z, Z_tv, t_tv, max_v, beta, update, update_rel, facility_idx, hosp_begin, max_d, parallel, nthreads);
    return R_NilValue;
END_RCPP
}
// update_alpha11
void update_alpha11(arma::colvec& alpha, arma::colvec& S0, arma::colvec& t_start, arma::colvec& t_end, arma::colvec& events_per_day, arma::mat& exp_z_beta0, arma::colvec& diffalpha, arma::colvec& alpha_star, arma::mat& z, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, arma::colvec& beta, arma::colvec& facility, int num_facility, int D, arma::colvec& lambda, int num_events, arma::colvec& events_per_facility, int N, arma::mat& Epre, arma::mat& hosp_begin, arma::colvec& max_d, bool parallel, const unsigned int& nthreads);
RcppExport SEXP _tdrecur_update_alpha11(SEXP alphaSEXP, SEXP S0SEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP events_per_daySEXP, SEXP exp_z_beta0SEXP, SEXP diffalphaSEXP, SEXP alpha_starSEXP, SEXP zSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP betaSEXP, SEXP facilitySEXP, SEXP num_facilitySEXP, SEXP DSEXP, SEXP lambdaSEXP, SEXP num_eventsSEXP, SEXP events_per_facilitySEXP, SEXP NSEXP, SEXP EpreSEXP, SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP parallelSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type events_per_day(events_per_daySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type exp_z_beta0(exp_z_beta0SEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type diffalpha(diffalphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha_star(alpha_starSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type num_events(num_eventsSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type events_per_facility(events_per_facilitySEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Epre(EpreSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nthreads(nthreadsSEXP);
    update_alpha11(alpha, S0, t_start, t_end, events_per_day, exp_z_beta0, diffalpha, alpha_star, z, Z_tv, t_tv, max_v, beta, facility, num_facility, D, lambda, num_events, events_per_facility, N, Epre, hosp_begin, max_d, parallel, nthreads);
    return R_NilValue;
END_RCPP
}
// cox_breslow_dep
List cox_breslow_dep(arma::mat& hosp_begin, arma::colvec& max_d, arma::colvec& t_start, arma::colvec& t_end, arma::mat& z, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, double tolb, double tola, int maxiter, arma::colvec& facility, int num_facility, int D, arma::colvec beta, bool parallel, const unsigned int& nthreads);
RcppExport SEXP _tdrecur_cox_breslow_dep(SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP zSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP tolbSEXP, SEXP tolaSEXP, SEXP maxiterSEXP, SEXP facilitySEXP, SEXP num_facilitySEXP, SEXP DSEXP, SEXP betaSEXP, SEXP parallelSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< double >::type tolb(tolbSEXP);
    Rcpp::traits::input_parameter< double >::type tola(tolaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(cox_breslow_dep(hosp_begin, max_d, t_start, t_end, z, Z_tv, t_tv, max_v, tolb, tola, maxiter, facility, num_facility, D, beta, parallel, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// compute_d12
List compute_d12(arma::mat& hosp_begin, arma::colvec& max_d, int D, arma::mat& z, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, arma::colvec& facility, int num_facility);
RcppExport SEXP _tdrecur_compute_d12(SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP DSEXP, SEXP zSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP facilitySEXP, SEXP num_facilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_d12(hosp_begin, max_d, D, z, Z_tv, t_tv, max_v, facility, num_facility));
    return rcpp_result_gen;
END_RCPP
}
// ddloglik_cpp6
void ddloglik_cpp6(double l, double lchange, arma::colvec& L1, arma::mat& L2, arma::colvec& t_start, arma::colvec& t_end, arma::rowvec& Sm, arma::mat events_per_day_facility, arma::mat& z, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, arma::colvec& beta, int N, List facility_idx, arma::mat& hosp_begin, arma::colvec& max_d);
RcppExport SEXP _tdrecur_ddloglik_cpp6(SEXP lSEXP, SEXP lchangeSEXP, SEXP L1SEXP, SEXP L2SEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP SmSEXP, SEXP events_per_day_facilitySEXP, SEXP zSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP betaSEXP, SEXP NSEXP, SEXP facility_idxSEXP, SEXP hosp_beginSEXP, SEXP max_dSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    Rcpp::traits::input_parameter< double >::type lchange(lchangeSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type L2(L2SEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type Sm(SmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type events_per_day_facility(events_per_day_facilitySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type facility_idx(facility_idxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type max_d(max_dSEXP);
    ddloglik_cpp6(l, lchange, L1, L2, t_start, t_end, Sm, events_per_day_facility, z, Z_tv, t_tv, max_v, beta, N, facility_idx, hosp_begin, max_d);
    return R_NilValue;
END_RCPP
}
// ddloglik_cpp14
void ddloglik_cpp14(arma::colvec& alpha, arma::colvec& t_start, arma::colvec& t_end, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, arma::mat& z, arma::colvec& beta, int D, int N, arma::colvec& S0, arma::mat& exp_z_beta0);
RcppExport SEXP _tdrecur_ddloglik_cpp14(SEXP alphaSEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP zSEXP, SEXP betaSEXP, SEXP DSEXP, SEXP NSEXP, SEXP S0SEXP, SEXP exp_z_beta0SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type exp_z_beta0(exp_z_beta0SEXP);
    ddloglik_cpp14(alpha, t_start, t_end, Z_tv, t_tv, max_v, z, beta, D, N, S0, exp_z_beta0);
    return R_NilValue;
END_RCPP
}
// compute_facility_idx_ind
List compute_facility_idx_ind(IntegerVector& facility, int num_facility);
RcppExport SEXP _tdrecur_compute_facility_idx_ind(SEXP facilitySEXP, SEXP num_facilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_facility_idx_ind(facility, num_facility));
    return rcpp_result_gen;
END_RCPP
}
// compute_d
List compute_d(NumericMatrix& hosp_begin, IntegerVector& max_d, int D, NumericMatrix& z, IntegerVector& facility, int num_facility);
RcppExport SEXP _tdrecur_compute_d(SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP DSEXP, SEXP zSEXP, SEXP facilitySEXP, SEXP num_facilitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_d(hosp_begin, max_d, D, z, facility, num_facility));
    return rcpp_result_gen;
END_RCPP
}
// ddloglik_cpp
void ddloglik_cpp(NumericVector& L1, NumericMatrix& L2, NumericVector& S0, IntegerVector& t_start, IntegerVector& t_end, NumericVector& Sm, IntegerMatrix events_per_day_facility, NumericVector& exp_z_beta, NumericMatrix& z, arma::colvec& beta, List facility_idx);
RcppExport SEXP _tdrecur_ddloglik_cpp(SEXP L1SEXP, SEXP L2SEXP, SEXP S0SEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP SmSEXP, SEXP events_per_day_facilitySEXP, SEXP exp_z_betaSEXP, SEXP zSEXP, SEXP betaSEXP, SEXP facility_idxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type L2(L2SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Sm(SmSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type events_per_day_facility(events_per_day_facilitySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type exp_z_beta(exp_z_betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type facility_idx(facility_idxSEXP);
    ddloglik_cpp(L1, L2, S0, t_start, t_end, Sm, events_per_day_facility, exp_z_beta, z, beta, facility_idx);
    return R_NilValue;
END_RCPP
}
// update_beta_ind
void update_beta_ind(NumericVector& L1, NumericMatrix& L2, NumericVector& S0, NumericVector& Sm, IntegerVector& t_start, IntegerVector& t_end, IntegerMatrix events_per_day_facility, NumericVector& exp_z_beta0, NumericMatrix& z, arma::colvec& beta, arma::colvec& update, arma::colvec& z_beta, IntegerVector& facility, int num_facility, int D, int N, List facility_idx);
RcppExport SEXP _tdrecur_update_beta_ind(SEXP L1SEXP, SEXP L2SEXP, SEXP S0SEXP, SEXP SmSEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP events_per_day_facilitySEXP, SEXP exp_z_beta0SEXP, SEXP zSEXP, SEXP betaSEXP, SEXP updateSEXP, SEXP z_betaSEXP, SEXP facilitySEXP, SEXP num_facilitySEXP, SEXP DSEXP, SEXP NSEXP, SEXP facility_idxSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type L2(L2SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Sm(SmSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type events_per_day_facility(events_per_day_facilitySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type exp_z_beta0(exp_z_beta0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type update(updateSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type z_beta(z_betaSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type facility_idx(facility_idxSEXP);
    update_beta_ind(L1, L2, S0, Sm, t_start, t_end, events_per_day_facility, exp_z_beta0, z, beta, update, z_beta, facility, num_facility, D, N, facility_idx);
    return R_NilValue;
END_RCPP
}
// ddloglik_Ui
arma::mat ddloglik_Ui(NumericVector& L1, NumericMatrix& L2, NumericVector& S0, IntegerVector& t_start, IntegerVector& t_end, NumericVector& Sm, IntegerMatrix events_per_day_facility, NumericVector& exp_z_beta, NumericMatrix& z, arma::colvec& beta, List facility_idx, NumericMatrix& hosp_begin, IntegerVector& max_d);
RcppExport SEXP _tdrecur_ddloglik_Ui(SEXP L1SEXP, SEXP L2SEXP, SEXP S0SEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP SmSEXP, SEXP events_per_day_facilitySEXP, SEXP exp_z_betaSEXP, SEXP zSEXP, SEXP betaSEXP, SEXP facility_idxSEXP, SEXP hosp_beginSEXP, SEXP max_dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type L2(L2SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type Sm(SmSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type events_per_day_facility(events_per_day_facilitySEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type exp_z_beta(exp_z_betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< List >::type facility_idx(facility_idxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type max_d(max_dSEXP);
    rcpp_result_gen = Rcpp::wrap(ddloglik_Ui(L1, L2, S0, t_start, t_end, Sm, events_per_day_facility, exp_z_beta, z, beta, facility_idx, hosp_begin, max_d));
    return rcpp_result_gen;
END_RCPP
}
// ddloglik_cpp11
void ddloglik_cpp11(NumericVector& S0, IntegerVector& t_start, IntegerVector& t_end, NumericVector& exp_z_beta);
RcppExport SEXP _tdrecur_ddloglik_cpp11(SEXP S0SEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP exp_z_betaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type exp_z_beta(exp_z_betaSEXP);
    ddloglik_cpp11(S0, t_start, t_end, exp_z_beta);
    return R_NilValue;
END_RCPP
}
// update_alpha_ind
void update_alpha_ind(NumericVector& alpha, NumericVector& S0, IntegerVector& t_start, IntegerVector& t_end, IntegerVector dm, NumericVector& exp_z_beta, NumericVector& exp_z_beta0, NumericVector& diffalpha, NumericVector& alpha_star, IntegerVector& facility, int num_facility, int D, NumericVector& lambda, int num_events, NumericVector& O, int N, NumericMatrix& Epre);
RcppExport SEXP _tdrecur_update_alpha_ind(SEXP alphaSEXP, SEXP S0SEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP dmSEXP, SEXP exp_z_betaSEXP, SEXP exp_z_beta0SEXP, SEXP diffalphaSEXP, SEXP alpha_starSEXP, SEXP facilitySEXP, SEXP num_facilitySEXP, SEXP DSEXP, SEXP lambdaSEXP, SEXP num_eventsSEXP, SEXP OSEXP, SEXP NSEXP, SEXP EpreSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type exp_z_beta(exp_z_betaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type exp_z_beta0(exp_z_beta0SEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type diffalpha(diffalphaSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type alpha_star(alpha_starSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type num_events(num_eventsSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type O(OSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type Epre(EpreSEXP);
    update_alpha_ind(alpha, S0, t_start, t_end, dm, exp_z_beta, exp_z_beta0, diffalpha, alpha_star, facility, num_facility, D, lambda, num_events, O, N, Epre);
    return R_NilValue;
END_RCPP
}
// cox_breslow_ind
List cox_breslow_ind(NumericMatrix& hosp_begin, IntegerVector& max_d, IntegerVector& t_start, IntegerVector& t_end, NumericMatrix& z, double tolb, double tola, int maxiter, IntegerVector& facility, int num_facility, int D, arma::colvec beta);
RcppExport SEXP _tdrecur_cox_breslow_ind(SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP zSEXP, SEXP tolbSEXP, SEXP tolaSEXP, SEXP maxiterSEXP, SEXP facilitySEXP, SEXP num_facilitySEXP, SEXP DSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< NumericMatrix& >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type tolb(tolbSEXP);
    Rcpp::traits::input_parameter< double >::type tola(tolaSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(cox_breslow_ind(hosp_begin, max_d, t_start, t_end, z, tolb, tola, maxiter, facility, num_facility, D, beta));
    return rcpp_result_gen;
END_RCPP
}
// compute_d_omp
List compute_d_omp(arma::mat& hosp_begin, arma::colvec& max_d, int D, arma::mat& z, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, arma::colvec& facility, int num_facility, const unsigned int& nthreads);
RcppExport SEXP _tdrecur_compute_d_omp(SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP DSEXP, SEXP zSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP facilitySEXP, SEXP num_facilitySEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type facility(facilitySEXP);
    Rcpp::traits::input_parameter< int >::type num_facility(num_facilitySEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_d_omp(hosp_begin, max_d, D, z, Z_tv, t_tv, max_v, facility, num_facility, nthreads));
    return rcpp_result_gen;
END_RCPP
}
// ddloglik_omp
void ddloglik_omp(arma::colvec& L1, arma::mat& L2, arma::colvec& t_start, arma::colvec& t_end, arma::rowvec& Sm, arma::mat events_per_day_facility, arma::mat& z, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, arma::colvec& beta, int N, List facility_idx, arma::mat& hosp_begin, arma::colvec& max_d, const unsigned int& nthreads);
RcppExport SEXP _tdrecur_ddloglik_omp(SEXP L1SEXP, SEXP L2SEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP SmSEXP, SEXP events_per_day_facilitySEXP, SEXP zSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP betaSEXP, SEXP NSEXP, SEXP facility_idxSEXP, SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type L1(L1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type L2(L2SEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< arma::rowvec& >::type Sm(SmSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type events_per_day_facility(events_per_day_facilitySEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type facility_idx(facility_idxSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nthreads(nthreadsSEXP);
    ddloglik_omp(L1, L2, t_start, t_end, Sm, events_per_day_facility, z, Z_tv, t_tv, max_v, beta, N, facility_idx, hosp_begin, max_d, nthreads);
    return R_NilValue;
END_RCPP
}
// ddloglik_all_omp
void ddloglik_all_omp(arma::colvec& alpha, arma::colvec& t_start, arma::colvec& t_end, arma::cube& Z_tv, arma::cube& t_tv, arma::mat& max_v, arma::mat& z, arma::colvec& beta, int D, int N, arma::colvec& S0, arma::mat& exp_z_beta0, arma::mat& hosp_begin, arma::colvec& max_d, const unsigned int& nthreads);
RcppExport SEXP _tdrecur_ddloglik_all_omp(SEXP alphaSEXP, SEXP t_startSEXP, SEXP t_endSEXP, SEXP Z_tvSEXP, SEXP t_tvSEXP, SEXP max_vSEXP, SEXP zSEXP, SEXP betaSEXP, SEXP DSEXP, SEXP NSEXP, SEXP S0SEXP, SEXP exp_z_beta0SEXP, SEXP hosp_beginSEXP, SEXP max_dSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_start(t_startSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type t_end(t_endSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type Z_tv(Z_tvSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type t_tv(t_tvSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type max_v(max_vSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type exp_z_beta0(exp_z_beta0SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type hosp_begin(hosp_beginSEXP);
    Rcpp::traits::input_parameter< arma::colvec& >::type max_d(max_dSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type nthreads(nthreadsSEXP);
    ddloglik_all_omp(alpha, t_start, t_end, Z_tv, t_tv, max_v, z, beta, D, N, S0, exp_z_beta0, hosp_begin, max_d, nthreads);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tdrecur_compute_facility_idx", (DL_FUNC) &_tdrecur_compute_facility_idx, 2},
    {"_tdrecur_update_beta", (DL_FUNC) &_tdrecur_update_beta, 20},
    {"_tdrecur_update_alpha11", (DL_FUNC) &_tdrecur_update_alpha11, 25},
    {"_tdrecur_cox_breslow_dep", (DL_FUNC) &_tdrecur_cox_breslow_dep, 17},
    {"_tdrecur_compute_d12", (DL_FUNC) &_tdrecur_compute_d12, 9},
    {"_tdrecur_ddloglik_cpp6", (DL_FUNC) &_tdrecur_ddloglik_cpp6, 17},
    {"_tdrecur_ddloglik_cpp14", (DL_FUNC) &_tdrecur_ddloglik_cpp14, 12},
    {"_tdrecur_compute_facility_idx_ind", (DL_FUNC) &_tdrecur_compute_facility_idx_ind, 2},
    {"_tdrecur_compute_d", (DL_FUNC) &_tdrecur_compute_d, 6},
    {"_tdrecur_ddloglik_cpp", (DL_FUNC) &_tdrecur_ddloglik_cpp, 11},
    {"_tdrecur_update_beta_ind", (DL_FUNC) &_tdrecur_update_beta_ind, 17},
    {"_tdrecur_ddloglik_Ui", (DL_FUNC) &_tdrecur_ddloglik_Ui, 13},
    {"_tdrecur_ddloglik_cpp11", (DL_FUNC) &_tdrecur_ddloglik_cpp11, 4},
    {"_tdrecur_update_alpha_ind", (DL_FUNC) &_tdrecur_update_alpha_ind, 17},
    {"_tdrecur_cox_breslow_ind", (DL_FUNC) &_tdrecur_cox_breslow_ind, 12},
    {"_tdrecur_compute_d_omp", (DL_FUNC) &_tdrecur_compute_d_omp, 10},
    {"_tdrecur_ddloglik_omp", (DL_FUNC) &_tdrecur_ddloglik_omp, 16},
    {"_tdrecur_ddloglik_all_omp", (DL_FUNC) &_tdrecur_ddloglik_all_omp, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_tdrecur(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
