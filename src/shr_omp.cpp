//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>
#include <omp.h>

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' created on 2/14/2022: this one use stratified model

//' this function calculates Sm (1 by p), sum of z in event set; events_per_facility, the vector of total events in each
//' facility; num_events, number of all events; events_per_day_facility, number of events each day in each facility
//'@param hosp_begin matrix of hospitalization begin date
//'@param max_d the largest number of hospitalization in each segment
//'@param D number of distinct failure time number of days)
//'@param z time-independent variables
// [[Rcpp::export]]
List compute_d_omp(arma::mat &hosp_begin, arma::colvec &max_d, int D,
arma::mat &z, arma::cube &Z_tv, arma::cube &t_tv, arma::mat &max_v,
	arma::colvec &facility, int num_facility,
                   const unsigned int &nthreads=1){
	  int N = hosp_begin.n_rows;
	  int p1 = z.n_cols;
    int p2 = Z_tv.n_slices;
    int p = p1 + p2;
    arma::rowvec Sm(p); //Sm vector
    arma::colvec events_per_facility(num_facility);
    int num_events = 0;
    arma::mat events_per_day_facility(num_facility, D);
    arma::colvec events_per_day(D);
    omp_set_num_threads(nthreads);
#pragma omp parallel
{
  arma::rowvec local_Sm(p);
  arma::colvec local_events_per_facility(num_facility);
  arma::mat local_events_per_day_facility(num_facility, D);
  int local_num_events = 0;
  arma::colvec local_events_per_day(D);
  arma::mat tempz, tempt;
#pragma omp for schedule(guided)
    // i for patients, s for parameters
    for(int i = 0; i < N; i ++){
    	if(max_d(i) > 0){
      	local_events_per_facility(facility(i)) += max_d(i);
      	local_num_events += max_d(i);
      	for (int s = 0; s < p1; s ++){
      		  local_Sm(s) += z(i, s) * max_d(i);
	  	  }

    	  for (int j = 0; j < max_d(i); j ++ ){
    	  	  local_events_per_day_facility(facility(i), hosp_begin(i, j) - 1) += 1;
    	  	  local_events_per_day(hosp_begin(i, j) - 1) += 1;
            for (int s = 0; s < p2; s++) {
              tempz = Z_tv.slice(s);
              tempt = t_tv.slice(s);

              for (int t = max_v(i, s) - 1; t >= 0; t--) {
                if (tempt(i, t) <= hosp_begin(i, j)) {
                  local_Sm(s + p1) += tempz(i, t + 1);
                  break;
                }
              }
              if (tempt(i, 0) > hosp_begin(i, j)) {
                local_Sm(s + p1) += tempz(i, 0);
              }
            }
    	  }
    	}
    }
    #pragma omp critical
{
  for(int d = 0; d < D; d++){
    events_per_day(d) += local_events_per_day(d);
    for(int f = 0; f < num_facility; f ++){
      events_per_day_facility(f, d) += local_events_per_day_facility(f, d);
    }
  }
  for(int s = 0; s < p; s ++){
    Sm(s) += local_Sm(s);
  }
  for(int f = 0; f < num_facility; f ++){
    events_per_facility(f) += local_events_per_facility(f);
  }
  num_events += local_num_events;

}
}
    return List::create(Named("Sm") = Sm,
    	Named("events_per_day_facility") = events_per_day_facility,
    	Named("events_per_day") = events_per_day,
    	Named("events_per_facility") = events_per_facility,
    	Named("num_events") = num_events);
}

//'@param l log partial likelihood
//'@param S0 D+1 vector, we will do day 0 (begining of day 1) to day D (end of day D)
//'@param t_start a vector of segment start date, starting from 1, minus 1 when using it.
//'@param t_end a vector of segment end date, starting from 1.
//'@param Sm 1 by p vector
//'@param ym number of changes on each day
// [[Rcpp::export]]
void ddloglik_omp(arma::colvec &L1, arma::mat &L2,
arma::colvec &t_start, arma::colvec &t_end, arma::rowvec &Sm,
arma::mat events_per_day_facility,
arma::mat& z, arma::cube &Z_tv, arma::cube &t_tv,
                          arma::mat &max_v,
arma::colvec &beta, int N, List facility_idx,
                   arma::mat &hosp_begin, arma::colvec &max_d,
                   const unsigned int &nthreads=1){
  L1.fill(0);
  L2.fill(0);
	//l = 0;
	//int N = t_start.n_rows;
	int D = events_per_day_facility.n_cols; //number of distinct failure time (number of days)
	int p1 = z.n_cols;
	int p2 = Z_tv.n_slices;
  int p = p1 + p2;
omp_set_num_threads(nthreads);
#pragma omp parallel
{
  arma::mat zt = arma::zeros<arma::mat>(N, p2);
  arma::colvec W;
  W.ones(N);
  arma::mat tempz, tempt;

  double S0 = 0.0;
  arma::rowvec S1(p);
  arma::mat S2(p, p);
  double exp_z_beta0;
  arma::colvec zvi_beta =
    z * beta.rows(0, p1 - 1); // zbeta for variables that are time-independent
  arma::rowvec zi(p);
  arma::mat zti_beta(1, 1); // zbeta for time t, patient i
#pragma omp for schedule(guided)
	// i for patients, j, s for covariates, k for time
  for (int k = 0; k < D; k++) { // for time

    zt.fill(0);
    for (int s = 0; s < p2; s++) { // for covariates
      tempz = Z_tv.slice(s); // now it's array, we can get one row
      tempt = t_tv.slice(s);
      for (int i = 0; i < N; i++) { // for patient
        // is it possible to combine these two for-loops and make it faster?
        // some steps below are not needed if all change times are between
        // t_start and t_end;
        // todo ???
        if ((W(i) == 1) & (k >= (t_start(i) - 1)) & (k < t_end(i))) {
          // only process when patient i is at risk on day k
          zt(i, s) = tempz(i, 0); // initialize the vector as 0s
          for (int flag = 0; flag < max_v(i, s); flag++) {
            if (tempt(i, flag) - 1 <= k) {
              // update if change happens before k
              zt(i, s) = tempz(i, flag + 1);
            }
          }
        }
      }
    }

  	for (int f=0; f<facility_idx.size(); f++)
  	{
      SEXP f_idx_temp = facility_idx[f];
      IntegerVector f_idx(f_idx_temp); //conversion
      int nj = f_idx.size();
  	  S0 = 0.0;
      S1.fill(0);
      S2.fill(0);

  		for(int i = 0; i < nj; i ++ ) {
    		if ((W(f_idx[i]) == 1) & (k >= (t_start(f_idx[i]) - 1)) & (k < t_end(f_idx[i]))) {
          zti_beta = zt.row(f_idx[i]) * beta.rows(p1, p - 1);
          exp_z_beta0 = exp(zvi_beta(f_idx[i]) + zti_beta(0, 0));
          S0 += exp_z_beta0;
          zi = join_rows(z.row(f_idx[i]), zt.row(f_idx[i]));
          S1 += zi * exp_z_beta0;
          S2 += zi.t() * zi * exp_z_beta0;
        }

    	}

      if(S0 == 0) // in case no one is at risk
      			{
      				S0 = 1e-6;
      			}
#pragma omp critical
{
      L2 += events_per_day_facility(f, k) / S0 * (S2 - S1.t() * S1 / S0);
      for (int s = 0; s < p; s++) {
        L1[s] -=
          events_per_day_facility(f, k) / S0 * S1(s);
      }

  	}
  }
}} // end parallel
	for(int j = 0; j < p1 + p2; j ++ ){
    	  L1[j] += Sm(j);
  }
}

//'This function calculates S0 for all patients
//'@param t_start a vector of segment start date, starting from 1, minus 1 when
// using it.
//'@param t_end a vector of segment end date, starting from 1.
//'@param t_tv list of time-dependent variables' changing time, all larger than
// 1, meaning not change on the first day
//'@param max_v max number of variable change, patient by variable
//'@param D number of distinct failure time (number of days)
//'@param exp_z_beta0 N by D matrix
// [[Rcpp::export]]
void ddloglik_all_omp(arma::colvec &alpha,
                   arma::colvec &t_start, arma::colvec &t_end,
                   arma::cube &Z_tv, arma::cube &t_tv,
                   arma::mat &max_v, arma::mat &z, arma::colvec &beta,
                   int D,
                   int N, arma::colvec &S0, arma::mat &exp_z_beta0,
                   arma::mat &hosp_begin, arma::colvec &max_d,
                   const unsigned int &nthreads=1) {

  // l = 0;
  int p1 = z.n_cols;
  int p2 = Z_tv.n_slices;
  int p = p1 + p2;

    omp_set_num_threads(nthreads);
#pragma omp parallel
{

  arma::mat zt = arma::zeros<arma::mat>(N, p2);
  arma::colvec W;
  W.ones(N);
  arma::mat tempz, tempt;
  double exp_z_beta;
  arma::colvec zvi_beta =
    z * beta.rows(0, p1 - 1); // zbeta for variables that are time-independent
  arma::mat zti_beta(1, 1); // zbeta for time t, patient i

  // i for patients, s for covariates, t for change point, j for time
#pragma omp for schedule(guided)
  for (int j = 0; j < D; j++) { // for time
    zt.fill(0);
    for (int s = 0; s < p2; s++) { // for covariates
      tempz = Z_tv.slice(s);
      tempt = t_tv.slice(s);
      for (int i = 0; i < N; i++) { // for patient
        // is it possible to combine these two for-loops and make it faster?
        // some steps below are not needed if all change times are between
        // t_start and t_end;
        // todo ???
        if ((W(i) == 1) & (j >= (t_start(i) - 1)) & (j < t_end(i))) {
          // only process when patient i is at risk on day j
          zt(i, s) = tempz(i, 0); // initialize the vector as 0s
          for (int flag = 0; flag < max_v(i, s); flag++) {
            if (tempt(i, flag) - 1 <= j) {
              // update if change happens before j
              zt(i, s) = tempz(i, flag + 1);
            }
          }
        }
      }
    }

    S0(j) = 0.0;
    for (int i = 0; i < N; i++) {
      if ((W(i) == 1) & (j >= (t_start(i) - 1)) & (j < t_end(i))) {
        zti_beta = zt.row(i) * beta.rows(p1, p - 1);
        exp_z_beta0(i, j) = exp(zvi_beta(i) + zti_beta(0, 0));
        exp_z_beta = exp(alpha(i)) * exp_z_beta0(i, j);
        S0(j) += exp_z_beta;
      }
    }

  }
}
}
