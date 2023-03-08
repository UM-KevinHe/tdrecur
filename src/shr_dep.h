#ifndef SHR_DEP_H
#define SHR_DEP_H

#include <Rcpp.h>
Rcpp::List compute_d12(arma::mat &hosp_begin, arma::colvec &max_d, int D,
arma::mat &z, arma::cube &Z_tv, arma::cube &t_tv, arma::mat &max_v,
	arma::colvec &facility, int num_facility);
void ddloglik_cpp6(double l, double lchange, arma::colvec &L1, arma::mat &L2,
arma::colvec &t_start, arma::colvec &t_end, arma::rowvec &Sm,
arma::mat events_per_day_facility,
arma::mat& z, arma::cube &Z_tv, arma::cube &t_tv,
                          arma::mat &max_v,
arma::colvec &beta, int N, List facility_idx,
                   arma::mat &hosp_begin, arma::colvec &max_d);
void ddloglik_cpp14(arma::colvec &alpha,
                   arma::colvec &t_start, arma::colvec &t_end,
                   arma::cube &Z_tv, arma::cube &t_tv,
                   arma::mat &max_v, arma::mat &z, arma::colvec &beta,
                   int D,
                   int N, arma::colvec &S0, arma::mat &exp_z_beta0);

#endif